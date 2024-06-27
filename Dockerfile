# build environment
FROM clux/muslrust:stable as builder

WORKDIR /home/rust

RUN cargo new --bin simmr

WORKDIR /home/rust/simmr

COPY simmr/ ./simmr
COPY simmrd/ ./simmrd
COPY shared/ ./shared
COPY Cargo.toml ./

RUN cargo build --release

# release layer
FROM alpine:latest

COPY --from=builder /home/rust/simmr/target/x86_64-unknown-linux-musl/release/simmr /usr/local/bin/simmr
#COPY --from=builder /home/rust/simmr/target/aarch64-unknown-linux-musl/release/simmr /usr/local/bin/simmr
#COPY --from=builder /home/rust/simmr/target/x86_64-unknown-linux-musl/release/simmr /usr/local/bin/simmr

CMD /bin/sh

