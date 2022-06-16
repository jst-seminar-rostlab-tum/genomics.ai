/* eslint-disable no-unused-vars,react/jsx-no-bind */
import React, { useCallback, useState, useEffect } from 'react';
import {
  Box,
  Grid,
  TextField,
  Button,
  Avatar,
  Snackbar,
  Alert,
  Typography,
} from '@mui/material';
import { useHistory, useLocation } from 'react-router-dom';
import NavBar from 'components/NavBar';
import Footer from 'components/Footer/old';
import styles from './passwordresetpage.module.css';
import logo from 'assets/logo.svg';
import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import Input from 'components/Input/Input';
import CustomButton from 'components/CustomButton/index';

function PasswordResetPage(props) {
  const [errors, setErrors] = useState({});
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);
  const [input, setInput] = useState({
    password: '',
    confirmpass: '',
  });

  const location = useLocation();

  const boxStyle = {
    position: 'relative',
    align: 'center',
    width: '60%',
    bgcolor: 'background.paper',
    boxShadow: 3,
    top: 32,
    p: 12,
    borderRadius: 3,
  };

  function validateInput() {
    if (input.password != null && input.confirmpass != null) {
      if (input.password === input.confirmpass) return true;
    }

    return false;
  }

  function doPWReset() {
    setErrors({});
    if (!validateInput()) {
      setErrors({ passwordconfirm: 'Password mismatch!' });
      return;
    }

    const params = new URLSearchParams(location.search);
    const token = params.get('token');
    // send POST request with token and password
    const requestOptions = {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        password: input.password,
      }),
    };
    fetch(`${BACKEND_ADDRESS}/password_reset/${token}`, requestOptions).then(
      (response) => {
        if (response.status === 404) {
          setErrors({ response: response.statusText });
        }
        if (response.status === 400) {
          setErrors({ response: response.statusText });
        }
      },
    );
    setSnackbarVisible(true);
    // history.push('/');
  }

  function onSnackbarClose() {
    setSnackbarVisible(false);
    setErrors({});
  }

  const handleTextChange = useCallback(
    (e) => {
      setInput((prevState) => ({
        ...prevState,
        [e.target.id]: e.target.value,
      }));
    },
    [setInput],
  );

  return (
    <div className={styles.headerContainer}>
      <NavBar />
      <Box justifyContent="center" display="flex">
        <Box sx={boxStyle}>
          <Grid>
            <Grid container direction="row" justifyContent="center">
              <Grid xs item />
              <Grid align="center">
                <Avatar src={logo} sx={{ width: 72, height: 72 }} />
                <h1>Password Reset</h1>
              </Grid>
              <Grid xs align="right" item />
            </Grid>
            <Input
              id="password"
              type="password"
              error={!!errors.password}
              helperText={errors.password}
              label="New Password"
              placeholder="Enter New Password"
              required
              onChangeEvent={handleTextChange}
            />
            <Input
              id="confirmpass"
              type="password"
              error={!!errors.passwordconfirm}
              helperText={errors.passwordconfirm}
              label="Confirm Password"
              placeholder="Reenter Password"
              required
              onChangeEvent={handleTextChange}
            />
            <Box mt={1}>
              <CustomButton
                type="primary"
                sx={{ ml: 1, width: '100%' }}
                onClick={doPWReset}
              >
                <Typography>Reset</Typography>
              </CustomButton>
            </Box>
          </Grid>
        </Box>
      </Box>
      <Snackbar
        open={isSnackbarVisible}
        autoHideDuration={3000}
        onClose={onSnackbarClose}
      >
        <Alert
          severity={errors.response ? 'error' : 'success'}
          sx={{ width: '100%' }}
          onClose={onSnackbarClose}
        >
          {errors.response ? errors.response : 'Password reset successfully!'}
        </Alert>
      </Snackbar>
    </div>
  );
}

export default PasswordResetPage;
