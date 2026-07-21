# Deployment and data-access guide

This guide covers local validation, container deployment, HTTPS, access control, monitoring, and safe atlas updates.

## 1. Before deployment

SpatialAtlas can be run on a workstation, an institutional server, or a cloud virtual machine. Before making an atlas reachable outside the local computer, confirm that the tissue images, expression values, coordinates, annotations, and citation text may be shown to the intended audience. The MIT licence applies to the software and does not grant permission to redistribute biological data.

Keep non-public data on loopback, a private subnet, or an institutional service protected by an appropriate identity layer. The browser receives expression vectors and spatial metadata for requested views, so hiding interface elements is not an access-control mechanism.

## 2. Prepare the deployment directory

Create a directory containing only the H5AD, image, marker, and configuration files used by the service. Check the production files before starting the container:

```bash
python validation/validate_artifacts.py \
  --h5ad /path/to/data/sample.h5ad \
  --image /path/to/data/tissue.png \
  --markers /path/to/data/markers.csv \
  --output-dir /path/to/validation \
  --prefix production
```

Give the service account read access but not write access to the biological files. Data should be mounted at runtime rather than copied into the application image.

## 3. Test on the local interface

```bash
cd software
docker compose config --quiet
docker build -t spatialatlas:production .
docker run --rm -p 127.0.0.1:8050:8050 \
  -v /absolute/path/to/data:/data:ro \
  -v /absolute/path/to/config.docker.yaml:/app/config.yaml:ro \
  spatialatlas:production
```

In a second terminal:

```bash
curl -fsS -o /dev/null -w '%{http_code}\n' http://127.0.0.1:8050/
curl -fsS -o /dev/null -w '%{http_code}\n' http://127.0.0.1:8050/_spatialatlas/tissue/0
```

Both requests should return 200. Also open the page locally and check one cluster view and one known spatial marker before configuring external access.

## 4. Docker Compose configuration

The supplied Compose file runs the application and an optional Caddy reverse proxy on the internal `web` network. Set deployment-specific paths as environment variables:

```bash
export SPATIALATLAS_DATA_DIR=/absolute/path/to/data
export SPATIALATLAS_CONFIG_FILE=/absolute/path/to/config.docker.yaml
export SPATIALATLAS_MEMORY_LIMIT=12g
cd software
docker compose config --quiet
```

Choose a memory limit above the measured viewer requirement with an allowance for concurrent requests. Data preparation, especially at 8 µm or finer resolution, is normally performed outside the application container; the compact H5AD is then mounted read-only for viewing.

## 5. HTTPS

Replace the placeholder domain in `Caddyfile` with a hostname whose DNS record points to the server. Permit inbound TCP 80 and 443 only where needed, then start the services:

```bash
docker compose up -d --build
docker compose ps
docker compose logs --tail=100 caddy spatialatlas
```

Caddy can request and renew a public certificate when DNS and ingress are configured correctly. Verify the proxy and application:

```bash
curl -fsS https://atlas.example.org/health
curl -fsS -o /dev/null -w '%{http_code}\n' https://atlas.example.org/
```

The supplied validation covers local container execution. Public HTTPS depends on the operator's domain and network configuration and therefore must be checked separately for each installation.

## 6. Access control and secrets

Publicly licensed data may be served without authentication if the data owner has approved this use. Restricted atlases should sit behind an institutional identity-aware proxy, VPN, private network, or equivalent authentication layer. Use a dedicated, unprivileged service account and read-only data mounts. Store credentials in a secrets manager or protected environment file; do not place passwords, tokens, private keys, or restricted download links in YAML files committed with the source.

## 7. Routine operation

Useful checks are:

```bash
docker compose ps
docker stats --no-stream
docker compose logs --since=1h spatialatlas caddy
curl -fsS https://atlas.example.org/health
```

The health endpoint confirms that the proxy can reach the application; it does not reopen every H5AD object or validate biological content. Include container restarts, HTTP errors, memory, disk space, and certificate renewal in the host's normal monitoring.

## 8. Updating an atlas

Keep versioned copies of the H5AD files, images, marker table, YAML configuration, checksum manifest, and container image tag. Prepare and check new files in a staging directory, run the viewer checks against a loopback port, and then update the read-only mount. Retain the preceding version until the replacement has been reviewed. If an atlas has been exposed unintentionally, remove external ingress first and follow the data owner's and institution's incident procedure.

## 9. Deployment checklist

| Item | Check |
|---|---|
| permission to display data | audience and reuse conditions confirmed |
| production files | artifact report completed and checksums stored |
| local service | atlas and tissue-image routes return HTTP 200 |
| file access | data and configuration mounted read-only |
| external service | valid HTTPS and suitable authentication/network restriction |
| operation | monitoring, backup, update, and rollback route documented |
