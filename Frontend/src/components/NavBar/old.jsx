import {
  Button, Toolbar, Typography, IconButton, Box,
} from '@mui/material';
import React, { useState, useCallback } from 'react';
import { Link } from 'react-router-dom';
import styles from './navbar.module.css';
import LoginForm from '../LoginForm';
import RegistrationForm from '../RegistrationForm';

import logo from '../../assets/logo.svg';

const NavBar = () => {
  // inserting the logic from the landing page here
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const onLoginClicked = useCallback(() => {
    setRegistrationFormVisible(false);
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

  const onSignUpClicked = useCallback(() => {
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible]);

  const onLoginFormClosed = useCallback(() => {
    setLoginFormVisible(false);
  }, [setLoginFormVisible]);

  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  // const { onLoginClicked, onSignUpClicked } = props;
  return (
    <div>
      <div className={styles.appbar}>
        <Toolbar>
          <Link to="/">
            <IconButton
              color="success"
              aria-label="open drawer"
            >
              <img alt="logo" src={logo} />
            </IconButton>
          </Link>
          <Typography sx={{ fontSize: '24px', fontWeight: '500' }}>genomics.ai</Typography>
          <Box sx={{ flexGrow: 1 }}>
            <Typography variant="h6" sx={{ flexGrow: 1 }}>
              <ul>
                <Box>
                  <li><Link to="/" className={styles.navbarcontent}>Home</Link></li>
                </Box>
                <li><Link to="about" className={styles.navbarcontent}>Team</Link></li>
                <li><Link to="docs" className={styles.navbarcontent}>Docs</Link></li>
                <li><Link to="contact" className={styles.navbarcontent}>Contact</Link></li>
              </ul>
            </Typography>
          </Box>
          <Box className={styles.buttonContainer}>
            <Button
              variant="text"
              color="inherit"
              sx={{ m: 2 }}
              onClick={onSignUpClicked}
            >
              <Typography sx={{ color: 'black', fontWeight: '500' }}>
                Sign up
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
        <LoginForm
          visible={isLoginFormVisible}
          onClose={onLoginFormClosed}
        />
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
