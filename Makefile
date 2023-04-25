.POSIX:

APP=fasta-api

test:
	func start

install:
	python3 -mvenv .venv
	. .venv/bin/activate \
	&& pip install --no-cache-dir -r requirements.txt

login:
	az login --use-device-code

publish:
	func azure functionapp publish $(APP)

