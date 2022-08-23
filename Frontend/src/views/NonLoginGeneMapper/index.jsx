import {
  useContext, React, useEffect, useState,
} from 'react';
import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import { useHistory } from 'react-router-dom';
import { Box, Snackbar, Alert } from '@mui/material';
import GeneMapper from 'views/GeneMapper';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import PasswordForgetForm from 'components/PasswordForgetForm';
import { LoginContext } from 'shared/context/loginContext';

// The non-login version of the genemapper page.
// The main genemapper component is used with loggedIn set to false.
function NonLoginGeneMapper() {
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  // get jwt token for temporary upload access
  useEffect(() => {
    // check if token exists && current token is not expired
    if (localStorage.getItem('jwt')) {
      const expiresAt = localStorage.getItem('temp_jwt_expiresAt');
      // check if the jwt token hasn't expired. Divide by 1000 since expiresAt is in seconds
      if (Math.ceil(Date.now() / 1000) < expiresAt) { return; }
    }

    fetch(`${BACKEND_ADDRESS}/temp_auth`)
      .then((response) => {
        if (!response.ok) {
          console.log('failed to get authorization for file upload.');
          setSnackbarVisible(true);
          return null;
        }
        return response.json();
      }) // set token
      .then((data) => {
        localStorage.setItem('jwt', data.jwt);
        // Add the timestamp when the jwt token was fetched.
        // Useful for checking if a token is expired
        localStorage.setItem('temp_jwt_expiresAt', data.expiresAt);
      });
  }, []);

  return (
    <Box sx={{
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'space-between',
      mineight: '100vh',
    }}
    >
      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          minHeight: '100vh',
        }}
      >
        {/* Non-login GeneMapper */}
        <GeneMapper sidebarShown={false} loggedIn={false} />
        {/* Snack bar for unsuccessful upload access */}
        <Snackbar
          open={isSnackbarVisible}
          autoHideDuration={10000}
          onClose={() => setSnackbarVisible(false)}
        >
          <Alert
            severity="error"
            sx={{ width: '100%' }}
            onClose={() => setSnackbarVisible(false)}
          >
            Failed to get upload authorization. Please retry later.
          </Alert>
        </Snackbar>
      </Box>
    </Box>
  );
}

export default NonLoginGeneMapper;
