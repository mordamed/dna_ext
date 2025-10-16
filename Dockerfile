# Simple demo for DNA Extension
FROM postgres:15

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    postgresql-server-dev-15 \
    && rm -rf /var/lib/apt/lists/*

# Copy extension source
COPY . /tmp/dna_ext/
WORKDIR /tmp/dna_ext

# Build and install extension
RUN make USE_PGXS=1 && make USE_PGXS=1 install

# Set up database
ENV POSTGRES_DB=dna_demo
ENV POSTGRES_USER=postgres
ENV POSTGRES_PASSWORD=password

# Copy demo script
COPY demo.sql /docker-entrypoint-initdb.d/