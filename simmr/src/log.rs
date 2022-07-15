/**
 * file: log.rs
 * desc: Application logging.
 */
use tracing::Level;
use tracing_subscriber;

/**
 * Sets up tracing and logging. By default all logging goes to stderr.
 */
pub fn setup_logging() {
    let _ = tracing_subscriber::fmt()
        .with_ansi(true)
        .with_target(false)
        .with_writer(std::io::stderr)
        .with_max_level(Level::INFO)
        .init();
}
