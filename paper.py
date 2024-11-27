import yaml
with open("pusar.yml". "r") as f:
	data = yaml.safe_load(f)
print(data)