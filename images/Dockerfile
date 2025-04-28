ARG VERSION=0.3.11

## BUILD PHASE #########################################################
FROM alpine:3.15 AS build-base

# Install build dependencies
RUN apk --no-cache add \
      build-base \
      bzip2-dev \
      curl-dev \
      libtbb-dev \
      ncurses-dev \
      python3-dev \
      wget \
      xz-dev \
      zlib-dev

# Build primer design program in virtualenv
FROM build-base AS build-primer-design
ARG VERSION

COPY crispr-primer-design /build

RUN python3 -m venv /opt/primer-design \
 && source /opt/primer-design/bin/activate \
 && pip3 install --upgrade pip setuptools \
 && pip3 install /build \
 && [[ "$(CRISPR_primer_design --version)" == "${VERSION}" ]]

# Download, patch and build Bowtie2 2.3.4.3
FROM build-base AS build-bowtie2

COPY <<'EOF' /tmp/bowtie2/bowtie2.patch
diff --color -ruN bowtie2-2.3.4.3/processor_support.h bowtie2-2.3.4.3.modified/processor_support.h
--- bowtie2-2.3.4.3/processor_support.h	2018-08-27 22:19:51.000000000 +0100
+++ bowtie2-2.3.4.3.modified/processor_support.h	2024-02-20 17:55:37.633253775 +0000
@@ -49,7 +49,7 @@
 #elif defined(USING_GCC_COMPILER)
         __get_cpuid(0x1, &regs.EAX, &regs.EBX, &regs.ECX, &regs.EDX);
 #else
-        std::cerr << “ERROR: please define __cpuid() for this build.\n”; 
+        std::cerr << "ERROR: please define __cpuid() for this build.\n"; 
         assert(0);
 #endif
         if( !( (regs.ECX & BIT(20)) && (regs.ECX & BIT(23)) ) ) return false;
EOF

RUN cd /tmp/bowtie2 \
 && wget -O bowtie2.zip https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.3/bowtie2-2.3.4.3-source.zip \
 && unzip bowtie2.zip \
 && patch -p0 < bowtie2.patch \
 && cd bowtie2-2.3.4.3 \
 && make

# Download, build and install Samtools 1.9
FROM build-base AS build-samtools

RUN mkdir -p /tmp/samtools \
 && cd /tmp/samtools \
 && wget -O samtools.tar.gz https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
 && tar -xf samtools.tar.gz \
 && cd samtools-1.9 \
 && ./configure --prefix=/opt/samtools \
 && make \
 && make install

# Download and install Bedtools 2.30.0
FROM build-base AS build-bedtools

RUN mkdir -p /opt/bedtools/bin \
 && cd /opt/bedtools/bin \
 && wget -O bedtools https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary \
 && chmod 755 bedtools

## RUNTIME #############################################################

FROM alpine:3.15
ARG VERSION

# Install runtime dependencies
# Nextflow expects Bash and GNU coreutils to exist
RUN apk add --no-cache \
      bash \
      coreutils \
      libbz2 \
      libcurl \
      libtbb \
      ncurses-libs \
      perl \
      python3 \
      xz-libs \
      zlib

# Install primer design program
COPY --from=build-primer-design /opt/primer-design /opt/primer-design

# Install Bowtie2
COPY --from=build-bowtie2 \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2 \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-align-s \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-align-l \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-build \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-build-s \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-build-l \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-inspect \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-inspect-s \
  /tmp/bowtie2/bowtie2-2.3.4.3/bowtie2-inspect-l \
  /opt/bowtie2/bin/

# Install Samtools
COPY --from=build-samtools /opt/samtools /opt/samtools

# Install Bedtools
COPY --from=build-bedtools /opt/bedtools /opt/bedtools

# Add binaries to PATH
ENV PATH $PATH:/opt/primer-design/bin:/opt/bowtie2/bin:/opt/samtools/bin:/opt/bedtools/bin

LABEL version="${VERSION}"
