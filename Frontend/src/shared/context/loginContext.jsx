import { createContext, useState } from "react";

const LoginContext = createContext()
LoginContext.displayName = "LoginContext"

const LoginProvider = (props) => {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);
  const [isPasswordForgetFormVisible, setPasswordForgetFormVisible] = useState(false);

  const onLoginFormClosed = () => {
    setLoginFormVisible(false);
  }

  const onRegistrationFormClosed = () => {
    setRegistrationFormVisible(false);
  }

  const onForgetPasswordClosed = () => {
    setPasswordForgetFormVisible(false);
  }

  const value = {
    loginVisible: isLoginFormVisible,
    registerVisible: isRegistrationFormVisible,
    forgetVisible: isPasswordForgetFormVisible,
    loginClose: onLoginFormClosed,
    registerClose: onRegistrationFormClosed,
    forgetClose: onForgetPasswordClosed,
    switchRegister: setRegistrationFormVisible,
    switchLogin: setLoginFormVisible,
    switchForget: setPasswordForgetFormVisible
  }

  return <LoginContext.Provider value={value} {...props}/>
}

export { LoginContext };
export default LoginProvider;