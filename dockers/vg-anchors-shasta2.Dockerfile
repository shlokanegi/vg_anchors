FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

# Minimal runtime libs for the one-file executable
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libstdc++6 \
    libgomp1 \
    libjansson4 \
    && rm -rf /var/lib/apt/lists/*

# Copy the prebuilt executable from the repository root build context
# Make sure you run docker build with the REPO ROOT as context
COPY ../dist/vg-anchors-0.1.0 /usr/local/bin/vg-anchors
RUN chmod +x /usr/local/bin/vg-anchors

# Default command
ENTRYPOINT ["/usr/local/bin/vg-anchors"]
