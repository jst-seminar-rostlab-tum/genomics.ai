// TODO: Delete this, we should not touch console object
export function init_environment() {
  let the_console: any = console;

  the_console.logCopy = console.log.bind(console);
  the_console.log = (...args: any[]) => {
    the_console.logCopy("[" + new Date().toUTCString() + "] ", ...args);
  };
  the_console.logWarnCopy = the_console.warn.bind(console);
  the_console.warn = (...args: any[]) => {
    the_console.logWarnCopy("[" + new Date().toUTCString() + "] ", ...args);
  };
  the_console.logErrCopy = the_console.error.bind(console);
  the_console.error = (...args: any) => {
    the_console.logErrCopy("[" + new Date().toUTCString() + "] ", ...args);
  };
}
