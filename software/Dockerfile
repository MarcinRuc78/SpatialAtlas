FROM python:3.11-slim

WORKDIR /app

COPY requirements-app.txt .
RUN pip install --no-cache-dir -r requirements-app.txt

COPY app.py .
COPY config.docker.yaml ./config.yaml
COPY assets ./assets

EXPOSE 8050
CMD ["gunicorn", "-b", "0.0.0.0:8050", "app:server", "--workers", "2", "--threads", "4", "--timeout", "120", "--preload"]
