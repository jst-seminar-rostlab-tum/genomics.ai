import { Box, Stack } from '@mui/material';
import React, { useContext, useEffect, useRef } from 'react';
import { useHistory } from 'react-router-dom';
import NavBar from 'components/NavBar';
import Header from 'components/general/Header';
import Footer from 'components/Footer';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import PasswordForgetForm from 'components/PasswordForgetForm';
import { LoginContext } from 'shared/context/loginContext';
import styles from './headerView.module.css';

// The GeneMapper HeaderView for both the logged-in and non-logged-in version.
// For the non-logged in version the navigation bar is rendered.
function HeaderView({
  title, rightOfTitle, replaceHeaderRight, children, loggedIn,
}) {
  const ref = useRef();

  useEffect(() => {
    console.log(ref);
  }, []);

  // The below variables are only necessary to initialize when
  // the non-logged in version is used.
  let onLoginClicked = null;
  let onSignUpClicked = null;
  let history = null;
  let context = null;

  // if not logged in, call the following code.
  if (!loggedIn) {
    context = useContext(LoginContext);

    onLoginClicked = () => {
      context.switchRegister(false);
      context.switchLogin(true);
    };

    onSignUpClicked = () => {
      context.switchLogin(false);
      context.switchRegister(true);
    };

    history = useHistory();
  }

  return (
    <Box
      ref={ref}
      sx={{
        width: '100%',
        height: '100vh',
        overflowY: 'auto',
      }}
    >
      {/* render navigation bar if not logged in */}
      {!loggedIn && (
        <NavBar
          onLoginClicked={onLoginClicked}
          onSignUpClicked={onSignUpClicked}
          executeScroll={() => history.push({ pathname: '/', state: { contact_us: true } })}
        />
      )}
      {/* render of necessary pop-ups when not logged in */}
      {context && context.loginVisible && <LoginForm />}
      {context && context.registerVisible && <RegistrationForm />}
      {context && context.forgetVisible && <PasswordForgetForm />}
      <Box sx={{
        minHeight: 'calc(100vh - 40px)',
        padding: '12px 60px 0 60px',
        boxSizing: 'border-box',
        display: 'flex',
      }}
      >
        <Stack
          direction="column"
          sx={{
            // height: "98vh",
            width: '100%',
            paddingLeft: '80px',
            flexGrow: 1,
          }}
        >
          <div style={{ paddingLeft: '60px', paddingRight: '60px' }}>
            <Header
              title={title}
              rightOfTitle={rightOfTitle}
              replaceRight={replaceHeaderRight}
              loggedIn={loggedIn}
            />
          </div>

          <Stack
            ref={ref}
            className="flexContainer"
            direction="row"
            sx={{ flexGrow: 1 }}
          >
            <div className={styles.content}>
              {children}
            </div>
          </Stack>
        </Stack>
      </Box>
      <Footer />
    </Box>
  );
}

export default HeaderView;
