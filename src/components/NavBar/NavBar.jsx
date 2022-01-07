import {
  AppBar, Button, Toolbar, Typography, IconButton, Box,
} from '@mui/material';
import React, { useState, useCallback } from 'react';
import styles from './navbar.module.css';
import LoginForm from '../LandingPage/LoginForm/LoginForm';
import RegistrationForm from '../LandingPage/RegistrationForm/RegistrationForm';

import logo from '../../assets/logo.svg';

const NavBar = () => {
  // inserting the logic from the landing page here
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

  // const { onLoginClicked, onSignUpClicked } = props;
  return (

    <div>
      <div className={styles.appbar}>
        <Toolbar>
          <IconButton color="success" aria-label="open drawer" href="/">
            <img src={logo} />
          </IconButton>
          <Typography sx={{ fontSize: '24px', fontWeight: '500' }}>genomics.ai</Typography>
          <Box sx={{ flexGrow: 1 }}>
            <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
              <ul>
                <Box>
                  <li><a href="/" className={styles.navbarcontent}>Home</a></li>
                </Box>
                <li><a href="about" className={styles.navbarcontent}>Team</a></li>
                <li><a href="docs" className={styles.navbarcontent}>Docs</a></li>
                <li><a href="contact" className={styles.navbarcontent}>Contact</a></li>
              </ul>
            </Typography>
          </Box>
          <Box className={styles.buttonContainer}>
            <Button variant="text" color="inherit" sx={{ m: 2 }} onClick={onSignUpClicked}>
              <Typography sx={{ color: 'black', fontWeight: '500' }}>
                Signup
              </Typography>
            </Button>
            <Button variant="outlined" color="info" onClick={onLoginClicked}>
              <Typography sx={{ color: '#4F83CC', fontWeight: '500' }}>
                Login
              </Typography>
            </Button>
          </Box>
        </Toolbar>
      </div>
      <div>
        <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />
        <RegistrationForm
          visible={isRegistrationFormVisible}
          onClose={onRegistrationFormClosed}
          onSuccessfulRegistration={onLoginClicked}
        />

      </div>
    </div>
  );
};

export default NavBar;

// <AppBar className={styles.appbar}> removed this from the above line before the beginning of the app
/*
<Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
          <ul>
            <Box>
              <li><a href="#section" className={styles.navbarcontent}>Home</a></li>
            </Box>
            <li><a href="#section" className={styles.navbarcontent}>Algorithm</a></li>
            <li><a href="#section" className={styles.navbarcontent}>References</a></li>
            <li><a href="#section" className={styles.navbarcontent}>Contact</a></li>
          </ul>
        </Typography>

*/
