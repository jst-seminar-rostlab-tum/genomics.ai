import {
  AppBar, Button, Toolbar, Typography, IconButton, Box,
} from '@mui/material';
import React, { useState, useCallback } from 'react';
import styles from './navbar.module.css';
import LoginForm from '../LandingPage/LoginForm/LoginForm';
import RegistrationForm from '../LandingPage/RegistrationForm/RegistrationForm';
// change logo to white
const Logo = () => (
  <svg width="56" height="56" viewBox="0 0 56 56" fill="white" xmlns="http://www.w3.org/2000/svg">
    <path
      d="M21.4739 39.0072C34.9113 41.2164 41.2622 35.2026 39.0055 21.4756L38.0237 17.095C37.3009 13.7017 36.6896 7.40488 39.5777 4.51687C40.358 3.73651 40.2955 2.40906 39.4387 1.55226C38.5819 0.695461 37.2545 0.632967 36.4747 1.41272C32.5961 5.29131 33.0049 12.1519 33.6592 16.1281L34.6273 20.5062C37.0609 32.1453 32.8829 35.6156 24.9001 35.2104C23.1504 35.1216 21.559 34.8489 20.505 34.6285L16.1276 33.6597C12.1508 33.006 5.29077 32.5967 1.41218 36.4753C0.632429 37.255 0.694925 38.5825 1.55172 39.4393C2.40852 40.2961 3.73597 40.3586 4.51572 39.5788C7.24607 36.8485 13.0458 37.19 17.0964 38.0371L21.4739 39.0072Z"
      fill="white"
    />
    <path
      d="M49.2983 21.1663C50.9424 20.6704 52.4738 19.8735 53.7052 18.6422C54.4849 17.8624 54.4224 16.535 53.5657 15.6782C52.7089 14.8214 51.3814 14.7589 50.6017 15.5387C50.0724 16.0679 49.4249 16.475 48.7045 16.7899L43.516 11.6014L40.0734 8.15876C39.6125 9.64841 39.5344 11.5368 39.6756 13.4099L40.6916 14.4259L43.9516 17.6858C42.6993 17.712 41.4468 17.6242 40.31 17.4786L41.1671 21.3191C41.1912 21.4642 41.208 21.6032 41.2301 21.7462C43.5253 21.965 46.2531 21.9661 48.7708 21.311C48.9547 21.2896 49.1302 21.239 49.2983 21.1663Z"
      fill="white"
    />
    <path
      d="M23.8989 20.7509L32.6381 29.4901C33.1038 28.1555 33.2376 26.4035 33.0298 24.2329L28.6116 19.8146C29.7295 19.8033 30.9586 19.8835 32.3354 20.0846L31.4974 16.2828C31.4732 16.1351 31.4502 15.9741 31.4269 15.8188C20.006 14.6617 14.6848 20.2179 15.8191 31.4267C15.9744 31.45 16.1361 31.4736 16.2831 31.4971L16.4064 31.5202L20.1088 32.3395C20.0219 31.677 19.9463 30.967 19.9077 30.2159C19.8808 29.6855 19.8761 29.1772 19.8847 28.6811L24.1754 32.9718C24.3855 32.9887 24.5905 33.0093 24.8061 33.0207C26.6924 33.1165 28.2157 32.9617 29.4219 32.5693L20.8049 23.9523C21.4469 22.5435 22.4482 21.4549 23.8989 20.7509Z"
      fill="white"
    />
    <path
      d="M21.1672 49.2947C21.2381 49.1285 21.2881 48.9549 21.311 48.7749C21.9667 46.2565 21.9655 43.5275 21.7467 41.231C21.6036 41.2089 21.4646 41.1921 21.3195 41.168L21.1962 41.1449L17.4686 40.32C17.6213 41.4981 17.7085 42.7378 17.6856 43.9518L13.4211 39.6873C11.5208 39.5383 9.63941 39.6161 8.15923 40.0743L16.7917 48.7067C16.475 49.4289 16.0664 50.0739 15.5391 50.6012C14.7593 51.381 14.8218 52.7084 15.6786 53.5652C16.5361 54.4227 17.8629 54.4845 18.6426 53.7048C19.874 52.4734 20.6713 50.9401 21.1672 49.2947Z"
      fill="white"
    />
  </svg>
);

const NavBar = (props) => {
  //inserting the logic from the landing page here
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

  //const { onLoginClicked, onSignUpClicked } = props;
  return (
    <div>
      <div className={styles.AppBar}>
      <Toolbar>
        <IconButton color="success" aria-label="open drawer">
          <a href="#section">
            <Logo className={styles.logo} />
          </a>
        </IconButton>
        <Typography className={styles.logoName}>GeneCruncher</Typography>
        <Box sx={{ flexGrow: 1 }}>
          <Typography variant="h6" component="div" sx={{ flexGrow: 1 }}>
            <ul>
              <Box>
                <li><a href="/" className={styles.navbarcontent}>Home</a></li>
              </Box>
              <li><a href="about" className={styles.navbarcontent}>About us</a></li>
              <li><a href="docs" className={styles.navbarcontent}>Docs</a></li>
              <li><a href="contact" className={styles.navbarcontent}>Contact</a></li>
            </ul>
          </Typography>
        </Box>

        <Button variant="text" color="inherit" sx={{ m: 2 }} onClick={onLoginClicked}>
          <Typography sx={{ color: 'black' }}>
            Login
          </Typography>
        </Button>
        <Button variant="outlined" color="info" onClick={onSignUpClicked}>Signup</Button>
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


//<AppBar className={styles.appbar}> removed this from the above line before the beginning of the app

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