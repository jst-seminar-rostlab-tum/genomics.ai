import NavBar from '../NavBar/NavBar'
import {useState} from "react";
import LoginForm from "./LoginForm/LoginForm";
import RegistrationForm from "./RegistrationForm/RegistrationForm";

function LandingPage() {

  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  function onLoginClicked() {
    setLoginFormVisible(true);
  }

  function onLoginFormClosed() {
    setLoginFormVisible(false);
  }

  function onSignUpClicked() {
    setRegistrationFormVisible(true);
  }

  function onRegistrationFormClosed() {
    setRegistrationFormVisible(false);
  }

  return (
    <div>
      <NavBar onLoginClicked={onLoginClicked} onSignUpClicked={onSignUpClicked}/>
      <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed}/>
      <RegistrationForm visible={isRegistrationFormVisible} onClose={onRegistrationFormClosed}
                        onSuccessfulRegistration={onLoginClicked}/>
    </div>
  );
}

export default LandingPage;