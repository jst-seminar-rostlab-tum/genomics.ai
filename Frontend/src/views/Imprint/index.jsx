import React, { useState, useCallback } from 'react';
import Box from '@mui/material/Box';
import Typography from '@mui/material/Typography';
import { useHistory } from 'react-router-dom';
import NavBar from 'components/NavBar';
import Footer from 'components/Footer';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';

export default function Imprint() {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const history = useHistory();

  const onLoginClicked = useCallback(() => {
    console.log('login');
    setRegistrationFormVisible(false);
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

  const onSignUpClicked = useCallback(() => {
    console.log('register');
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible]);

  const onLoginFormClosed = useCallback(() => {
    setLoginFormVisible(false);
  }, [setLoginFormVisible]);

  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  const executeScroll = () => history.push({ pathname: '/', state: { contact_us: true } });

  const regForm = isRegistrationFormVisible
    && <RegistrationForm visible={isRegistrationFormVisible} onClose={onRegistrationFormClosed} />;
  return (
    <Box sx={{ height: '100vh' }}>
      {isLoginFormVisible && <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />}
      {regForm}

      <Box>
        <NavBar
          position="relative"
          onLoginClicked={onLoginClicked}
          onSignUpClicked={onSignUpClicked}
          executeScroll={executeScroll}
        />
      </Box>

      <Box sx={{
        margin: '3em auto',
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        width: {
          xs: '90%', sm: '90%', md: '70%', lg: '70%', xl: '70%',
        },
      }}
      >
        <Typography fontWeight="bold" fontSize="1.4em">Imprint</Typography>
      </Box>

      <Box sx={{
        margin: '3em auto',
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        width: {
          xs: '90%',
          sm: '90%',
          md: '70%',
          lg: '70%',
          xl: '70%',
        },
      }}
      >
        <Typography fontWeight="bold" fontSize="1em">
          Vereinssitz und Geschäftsstelle Bonn
          <br />
          Ahrstraße 45, 53175 Bonn
          <br />
          Telefon: +49 (0) 228 30818-0
          <br />
          Telefax: +49 (0) 228 30818-30
          <br />
          Email: info (at) helmholtz.de
          <br />
          Webadresse: www.helmholtz.de
        </Typography>
      </Box>

      <div style={{ height: 'calc(100vh - 485px)' }} />

      <Footer />
    </Box>
  );
}
