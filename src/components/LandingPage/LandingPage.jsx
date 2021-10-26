import React, { useCallback, useState } from 'react';
import NavBar from '../NavBar/NavBar';
import LoginForm from './LoginForm/LoginForm';
import RegistrationForm from './RegistrationForm/RegistrationForm';

function LandingPage() {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const onLoginClicked = useCallback(() => {
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

  const onLoginFormClosed = useCallback(() => {
    setLoginFormVisible(false);
  }, [setLoginFormVisible]);

  const onSignUpClicked = useCallback(() => {
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible]);

  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  return (
    <div>
      <NavBar onLoginClicked={onLoginClicked} onSignUpClicked={onSignUpClicked} />
      <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />
      <RegistrationForm
        visible={isRegistrationFormVisible}
        onClose={onRegistrationFormClosed}
        onSuccessfulRegistration={onLoginClicked}
      />
    </div>
  );
}

export default LandingPage;
