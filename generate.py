#!/usr/bin/env python3
"""Filing-aware Walker constellation and optical ISL contact-plan generator.
Two-body circular orbital propagation. All generated ODC geometry is synthetic.
"""
import argparse, csv, json, math
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

RE=6378.137  # km, spherical Earth model
MU=398600.4418 # km^3/s^2
C=299792.458 # km/s
REGISTRY='https://github.com/irtf-spacerg/id-leo-constellations'


def write_csv(path,fields,rows):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def construct(n,planes,h,inc,phase,walker):
    if n%planes: raise ValueError('satellites must be divisible by planes')
    if not (0 <= inc <= 180): raise ValueError('inclination must be in [0,180]')
    s=n//planes
    if phase<0 or phase>=planes: raise ValueError('phase must be between 0 and planes-1')
    nodes=[]
    for p in range(planes):
        raan=360*p/planes*(0.5 if walker=='S' else 1.0)
        for k in range(s):
            u=(360*k/s+360*phase*p/n)%360
            nodes.append({'sat_id':f'SAT-{p:03d}-{k:03d}','plane':p,'slot':k,'altitude_km':h,'inclination_deg':inc,'raan_deg':raan,'argument_latitude_deg':u,'eccentricity':0.0})
    return nodes


def positions(nodes,t):
    h=np.array([v['altitude_km'] for v in nodes],float)
    inc=np.deg2rad([v['inclination_deg'] for v in nodes]); ra=np.deg2rad([v['raan_deg'] for v in nodes])
    u=np.deg2rad([v['argument_latitude_deg'] for v in nodes]); r=RE+h
    
    'For a 550 km circular orbit, the orbital period is approximately 95.6 minutes.'
    u=u+np.sqrt(MU/r**3)*t

    '''To convert to Earth-centered inertial coordinates (ECI), we apply two rotations: 
        i: rotates by the inclination. 
        O: rotates by the RAAN. 
        The final coordinates are given by the following equations:'''
    cu,su=np.cos(u),np.sin(u); 
    cO,sO=np.cos(ra),np.sin(ra);
    ci,si=np.cos(inc),np.sin(inc)

    'The returned array has shape: (len(nodes), 3), where each row corresponds to the (x, y, z) coordinates of a satellite in ECI (N,3)'
    return np.column_stack((r*(cO*cu-sO*su*ci),r*(sO*cu+cO*su*ci),r*su*si))


def feasible_edges(x,max_range):
    'We want to identify satellite pairs that could establish an optical inter-satellite link (ISL).'
    # Earth-clearance for each segment: min distance from origin to closed segment.
    # Vectorized by source to limit memory; still O(N^2) computational work.
    n=len(x)
    for i in range(n):
        d=x[i+1:]-x[i]
        distance=np.linalg.norm(d,axis=1) # The Euclidean distance between satellite i and all satellites j>i. axis=1 means calculate one distance for each row.
        valid=(distance>0)&(distance<=max_range) # Keeps only positive distances within the maximum assumed ISL range.
        if not valid.any():continue
        
        jidx=np.flatnonzero(valid)+i+1
        delta=x[jidx]-x[i]
        # calculates the closest point on the line segment joining the two satellites to Earth center.'
        lam=np.clip(-np.sum(x[i]*delta,axis=1)/np.sum(delta*delta,axis=1),0,1)

        # calculates the closest points to Earth along each candidate segment.'
        close=x[i]+lam[:,None]*delta

        # checks whether those closest points are outside Earth's spherical boundary.
        clear=np.linalg.norm(close,axis=1)>RE
        for j,dist in zip(jidx[clear],distance[jidx[clear]-i-1]):
            # yields the satellite pair (i,j) and their distance if they can establish a link without Earth obstruction.
            # yield instead of return allows this function to be used as an iterator, producing one result at a time.
            yield i,int(j),float(dist)


def schedule(n,edges,terminals,previous,hold_steps,step_idx,acq_steps):
    # Greedy nearest-neighbor selection; respects terminal degree. Links are
    # undirected optical connections; hold and acquisition are simplified.
    degree=[0]*n;selected={};used=set()

    '''
        sorted() orders the feasible edges according to a two-part key.
        First, links that were selected during the previous simulation step are prioritized.
        Second, links are sorted by distance.

        This is because the scheduler prioritizes existing connections over introducing new ones, even if a new link is shorter. 
        This approach helps maintain stable connections and reduces the frequency of link changes, which is important for optical inter-satellite links
          that require precise pointing and alignment.
    '''
    order=sorted(edges,key=lambda x:(0 if (x[0],x[1]) in previous else 1,x[2]))
    for i,j,d in order:
        key=(i,j)
        if degree[i]>=terminals or degree[j]>=terminals:continue
        if key not in previous and (i in used or j in used):continue
        selected[key]=d;degree[i]+=1;degree[j]+=1
    # Keep incumbent edges during minimum hold, displacing new selections if needed.
    mapping={(i,j):d for i,j,d in edges}
    for key,start in previous.items():
        if step_idx-start>=hold_steps or key not in mapping or key in selected:continue
        i,j=key
        for node in (i,j):
            if degree[node]>=terminals:
                victims=[k for k in selected if node in k and k not in previous]
                if victims:
                    victim=victims[-1]; del selected[victim]
                    for a in victim:degree[a]-=1
        if degree[i]<terminals and degree[j]<terminals:
            selected[key]=mapping[key];degree[i]+=1;degree[j]+=1

            """
            The result is a greedy degree-constrained graph. 
            It does not guarantee that the overall selected network is connected, that aggregate capacity is maximized, 
            or that AI workloads can communicate end to end.
            """
    return selected
"""
However, this implementation is only best-effort. If a previously selected link becomes geometrically infeasible, it cannot be maintained. 
Moreover, the scheduler first prioritizes all existing links, so much of the holding-time logic becomes redundant.
"""


def run(a):
    """
    - For starlink-gen1, the code uses a predefined shell obtained from the registry, including literature-assumed phasing.
    - For spacex-odc, it uses your supplied orbital parameters, checks whether the altitude lies between 500 and 2,000 km, and checks the inclination setting.
    """
    if a.preset=='starlink-gen1':
        n,p,h,inc,phase,walker=1584,72,550.0,53.0,39,'D'
        provenance='Starlink Gen1 shell gen1-s1; phasing=39 is registry-documented literature assumption'
    else:
        n,p,h,inc,phase,walker=a.satellites,a.planes,a.altitude,a.inclination,a.phasing,a.walker
        if not 500<=h<=2000:raise ValueError('ODC altitude must be between 500 and 2000 km')
        if abs(inc-30)>1e-9 and not a.assume_sso:
            raise ValueError('ODC requires inclination 30 or --assume-sso for an explicitly assumed sun-synchronous configuration')
        provenance='Synthetic Walker slice constrained by ODC filing altitude and inclination envelope; not an actual SpaceX design'
    if not (a.dt>0 and a.duration>0 and a.terminals>0 and a.rate_gbps>0 and a.max_range_km>0):raise ValueError('all time, terminal, rate and range parameters must be positive')
    nodes=construct(n,p,h,inc,phase,walker)

    if a.propagator=='sgp4':
        from orbit_sgp4 import (
            build_sgp4_satellites,
            sgp4_positions
        )

        epoch = datetime(
            2026, 10, 22,
            tzinfo=timezone.utc
        )

        satellites = build_sgp4_satellites(nodes, epoch)

    out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'satellites.csv',list(nodes[0]),nodes)

    # Geometry files: time-series positions; candidate/scheduled contacts as interval records.
    steps=int(math.ceil(a.duration/a.dt))
    pos_fields=['time_s','sat_id','x_km','y_km','z_km']
    link_fields=['src','dst','t_start_s','t_end_s','rate_bps','prop_delay_ms','distance_km']
    positions_writer=open(out/'positions.csv','w',newline='')
    pw=csv.DictWriter(positions_writer,fieldnames=pos_fields);pw.writeheader()

    # Maintain contact segments across sampling instants; closing at the next sampling point.
    active_candidate={};active_schedule={};previous={};
    candidate_rows=[];scheduled_rows=[]
    
    

    def update(active, now, current, rows,
            track_distance=True):

        def close_interval(k, old):
            rows.append({
                'src': nodes[k[0]]['sat_id'],
                'dst': nodes[k[1]]['sat_id'],
                't_start_s': old[0],
                't_end_s': now,
                'rate_bps': int(a.rate_gbps * 1e9),
                'prop_delay_ms': round(
                    old[1] / C * 1000, 5
                ),
                'distance_km': round(old[1], 4)
            })

        # Close links that disappeared or whose
        # propagation delay needs updating.
        for k, old in list(active.items()):
            if k not in current:
                close_interval(k, old)
                del active[k]

            elif track_distance:
                new_distance = current[k]
                old_distance = old[1]

                if abs(new_distance - old_distance) > 1e-6:
                    close_interval(k, old)
                    active[k] = (now, new_distance)

        # Start newly active links.
        for k, d in current.items():
            if k not in active:
                active[k] = (now, d)

    
    for step in range(steps):
        t=step*a.dt; 
        if a.propagator=='sgp4':
            xyz = sgp4_positions(satellites, epoch, t)
        else:
            xyz=positions(nodes,t)
        for nd,pos in zip(nodes,xyz):pw.writerow({'time_s':round(t,6),'sat_id':nd['sat_id'],'x_km':round(pos[0],5),'y_km':round(pos[1],5),'z_km':round(pos[2],5)})
        edges=list(feasible_edges(xyz,a.max_range_km))
        all_contacts={(i,j):d for i,j,d in edges}
        update(active_candidate,t,all_contacts,candidate_rows)
        selected=schedule(n,edges,a.terminals,previous,max(1,math.ceil(a.min_hold_s/a.dt)),step,0)
        
        # Simulated reacquisition: newly selected links do not transmit for this number of sampled steps.
        for key in selected:
            if key not in previous:previous[key]=step
        previous={k:previous[k] for k in selected}
        operational={k:d for k,d in selected.items() if step-previous[k]>=math.ceil(a.acquisition_s/a.dt)}
        update(active_schedule,t,operational,scheduled_rows, track_distance=True)
        if a.verbose:print(f't={t:g}s candidates={len(edges)} scheduled={len(selected)} ready={len(operational)}')
    
    
    positions_writer.close()
    update(active_candidate,a.duration,{},candidate_rows)
    update(active_schedule,a.duration,{},scheduled_rows)
    write_csv(out/'candidate_contacts.csv',link_fields,candidate_rows)
    write_csv(out/'scheduled_contacts.csv',link_fields,scheduled_rows)
    # Canonical directed trace includes both forward/backward directions and outage intervals.
    # Linux networking and tc/netem operate on directed traffic flows. Therefore, the code creates both directions for every scheduled link.
    trace_fields=['scenario_id','src','dst','t_start_s','t_end_s','available','rate_bps','prop_delay_ms']
    by_edge={}
    for row in scheduled_rows:
        for src,dst in [(row['src'],row['dst']),(row['dst'],row['src'])]:
            by_edge.setdefault((src,dst),[]).append(row)
    trace=[]
    for (src,dst),intervals in sorted(by_edge.items()):
        cursor=0
        for e in sorted(intervals,key=lambda r:r['t_start_s']):
            beg=e['t_start_s'];end=e['t_end_s']
            if beg>cursor:trace.append({'scenario_id':a.preset,'src':src,'dst':dst,'t_start_s':cursor,'t_end_s':beg,'available':0,'rate_bps':0,'prop_delay_ms':0})
            trace.append({'scenario_id':a.preset,'src':src,'dst':dst,'t_start_s':beg,'t_end_s':end,'available':1,'rate_bps':e['rate_bps'],'prop_delay_ms':e['prop_delay_ms']})
            cursor=end
        if cursor<a.duration:trace.append({'scenario_id':a.preset,'src':src,'dst':dst,'t_start_s':cursor,'t_end_s':a.duration,'available':0,'rate_bps':0,'prop_delay_ms':0})
    # Also, the file includes only directed links that were scheduled at least once. It does not create outage-only rows for every possible satellite pair.
    write_csv(out/'directed_trace.csv',trace_fields,trace)
    manifest={'scenario':a.preset,'generated_at_utc':datetime.now(timezone.utc).isoformat(),'registry_url':REGISTRY,'provenance':provenance,'registry_record_sha':{'starlink-gen1':'8fbf37e48f309d51f67f91d3bc6db234fa486adb','spacex-odc':'b7561294fb870f5b8729e80b9a1e5e52a96851e8'}[a.preset], 'orbital_model':'circular two-body ECI, fixed RAAN, no J2 drift, Earth treated as sphere for LOS','parameters':{'satellites':n,'planes':p,'satellites_per_plane':n//p,'altitude_km':h,'inclination_deg':inc,'walker':walker,'phase':phase,'phase_assumed':True,'duration_s':a.duration,'step_s':a.dt,'max_range_km':a.max_range_km,'terminals_per_satellite':a.terminals,'acquisition_s':a.acquisition_s,'min_hold_s':a.min_hold_s,'rate_gbps':a.rate_gbps},'limitations':['All optical link parameters and scheduling policies are scenario assumptions, not SpaceX specifications.','Schedules are sampled approximations, not exact continuous-time optical contact intervals.','Link propagation delay is approximated using distance at contact start; actual delay changes within a contact.','Minimum hold is best-effort under a greedy scheduler, and may break if geometry becomes infeasible.','Sun-synchronous orbits are not faithfully modeled by the two-body propagator.','No pointing angles, laser link budget, J2 perturbation or physical acquisition model.','Trace covers links seen as scheduled at least once, not all possible pairs.']}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    print(f'Wrote {out}: {n} satellites, {len(candidate_rows)} candidate intervals, {len(scheduled_rows)} scheduled intervals, {len(trace)} directed trace intervals')

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--preset',choices=['spacex-odc','starlink-gen1'],default='spacex-odc')
    ap.add_argument('--satellites',type=int,default=60)
    ap.add_argument('--planes',type=int,default=6)
    ap.add_argument('--altitude',type=float,default=550)
    ap.add_argument('--inclination',type=float,default=30)
    ap.add_argument('--phasing',type=int,default=1)
    ap.add_argument('--walker',choices=['D','S'],default='D')
    ap.add_argument('--assume-sso',action='store_true')
    ap.add_argument('--duration',type=float,default=1800,help='simulation seconds')
    ap.add_argument('--dt',type=float,default=30,help='sampling seconds')
    ap.add_argument('--terminals',type=int,default=4)
    ap.add_argument('--min-hold-s',type=float,default=0)
    ap.add_argument('--acquisition-s',type=float,default=0)
    ap.add_argument('--max-range-km',type=float,default=3500)
    ap.add_argument('--rate-gbps',type=float,default=10)
    ap.add_argument('--output',default='output')
    ap.add_argument('--verbose',action='store_true')
    ap.add_argument('--propagator', choices=['circular','sgp4'],default='circular',help='orbital two-body propagator; SGP4 is more accurate but slower')
    run(ap.parse_args())
