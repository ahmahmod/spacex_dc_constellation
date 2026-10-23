python3 -m venv .venv
source .venv/bin/activate

# pip install --upgrade pip
pip install -r requirements.txt

python generate.py \
  --preset spacex-odc \
  --satellites 600 \
  --planes 20 \
  --altitude 550 \
  --inclination 30 \
  --phasing 1 \
  --duration 3600 \
  --dt 30 \
  --terminals 4 \
  --output output/circular \
  --propagator circular

rm -rf .venv