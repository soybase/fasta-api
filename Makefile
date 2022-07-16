.POSIX:

run:
	source venv/bin/activate && uvicorn main:app

install:
	python3 -mvenv venv
	source venv/bin/activate \
	&& pip install --upgrade --no-cache-dir pip \
	&& pip install --no-cache-dir wheel \
	&& pip install --no-cache-dir -r requirements.txt