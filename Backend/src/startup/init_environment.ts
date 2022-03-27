export function init_environment(){
    let the_console : any = console;

    the_console.logCopy = console.log.bind(console);
    the_console.log = (data : any)=>{the_console.logCopy('[' + new Date().toUTCString() + '] ', data)};
    the_console.logWarnCopy = the_console.warn.bind(console);
    the_console.warn = (data : any)=>{the_console.logWarnCopy('[' + new Date().toUTCString() + '] ', data)};
    the_console.logErrCopy = the_console.error.bind(console);
    the_console.error = (data : any)=>{the_console.logErrCopy('[' + new Date().toUTCString() + '] ', data)};
}