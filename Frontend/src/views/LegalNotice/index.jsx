/* eslint-disable max-len */
import React, { useState, useCallback } from 'react';
import Box from '@mui/material/Box';
import Link from '@mui/material/Link';
import Typography from '@mui/material/Typography';
import { useHistory } from 'react-router-dom';
import NavBar from 'components/NavBar';
import Footer from 'components/Footer';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import { useAuth } from 'shared/context/authContext';

export default function LegalNotice() {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const history = useHistory();
  const [user, setUser] = useAuth();

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

  const executeScroll = () => user ? history.push({ pathname: '/sequencer/help'}) : history.push({ pathname: '/', state: { contact_us: true } });

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
        <Typography fontWeight="bold" fontSize="1.4em">Legal Notice</Typography>
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
        <Typography fontSize="1em">
          Helmholtz Zentrum München
          <br />
          Deutsches Forschungszentrum für Gesundheit und Umwelt (GmbH)
          <br />
          Ingolstädter Landstraße 1,
          <br />
          D-85764 Neuherberg, Germany
          <br />
          <br />
          Phone: +49 (0) 89 3187-0
          <br />
          email: info(at)helmholtz-muenchen.de
          <br />
          <br />
          <b>Representatives</b>
          <br />
          Board of Representatives: Prof. Dr. Matthias H. Tschöp, Kerstin Günther, Daniela Sommer (acting)
          <br />
          <br />
          <b>Chairwoman of Supervisory Board:</b>
          <br />
          MinDir’in Prof. Dr. Veronika von Messling
          <br />
          <br />
          <b>Register of Societies:</b>
          <br />
          Amtsgericht München HRB 6466
          <br />
          <br />
          <b>VAT ID number in accordance with § 27 a Umsatzsteuergesetz (German Turnover-Tax Law):</b> 
          <br />
          DE 129521671
          <br />
          <br />
          <br />
        </Typography>
      </Box>

      <div style={{ height: 'calc(100vh - 485px)' }} />

      <Footer />
    </Box>
  );
}
