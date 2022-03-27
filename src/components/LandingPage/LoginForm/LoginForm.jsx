/* eslint-disable */
import {
  Alert,
  Avatar, Box, Checkbox, FormControlLabel, Grid, Modal, Snackbar, TextField, Typography, Link,
} from '@mui/material';
import LoadingButton from '@mui/lab/LoadingButton';
import React, { useCallback, useState } from 'react';
import validator from 'validator';
import CloseIcon from '@mui/icons-material/Close';
import { useHistory } from 'react-router-dom';
import logo from '../../../assets/logo.svg';
import styles from './loginform.module.css';
import { BACKEND_ADDRESS } from '../../common/constants';
import PasswordForgetForm from '../PasswordReset/PasswordForgetForm';

function LoginForm(props) {
  const { setUser } = props;
  const [loginDetails, setLoginDetails] = useState({
    email: '',
    password: '',
    remember: false,
  });

  const [forgotPassword, setForgotPassword] = useState(false);
  const [errors, setErrors] = useState({});
  const [loading, setLoading] = useState(false);
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  const history = useHistory();

  const onClose = useCallback(() => {
    setLoginDetails({
      email: '',
      password: '',
      remember: false,
    });
    setErrors({});
    setLoading(false);
    props.onClose();
  }, [setLoginDetails, setErrors, setLoading, props]);

  const handleTextChange = useCallback((e) => {
    setLoginDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [setLoginDetails]);

  const handleCheckedChange = useCallback((e) => {
    setLoginDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.checked }));
  }, [setLoginDetails]);

  const handlepasswordForget = useCallback(() => {
    onClose();
    setForgotPassword(true);
  }, [setForgotPassword]);

  function validateInput() {
    let currentErrors = {};
    if (!validator.isEmail(loginDetails.email)) {
      currentErrors = { ...currentErrors, email: 'A valid e-mail is required!' };
    }
    if (loginDetails.password === '') {
      currentErrors = { ...currentErrors, password: 'Password is required!' };
    }
    setErrors(currentErrors);
    return !Object.keys(currentErrors).length;
  }

  async function onSuccessfulLogin(data) {
    localStorage.setItem('jwt', data.jwt);
    localStorage.setItem('user', JSON.stringify(data.user));
    onClose();
    await setUser(data.user);
    history.push('/sequencer/dashboard');
  }

  function onFailedLogin(code) {
    switch (code) {
      case 401:
        setErrors((prevState) => ({ ...prevState, response: 'Wrong credentials!' }));
        break;
      case 404:
        setErrors((prevState) => ({ ...prevState, response: 'User not found!' }));
        break;
      default:
        setErrors((prevState) => ({ ...prevState, response: 'Unknown error, please try again later!' }));
        break;
    }
  }

  const doLogin = useCallback(() => {
    if (!validateInput()) {
      return;
    }
    setLoading(true);
    // server communication
    const loginRequest = {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        Authorization: 'Bearer',
      },
      body: JSON.stringify({
        email: loginDetails.email,
        password: loginDetails.password,
      }),
    };
    fetch(`${BACKEND_ADDRESS}/auth`, loginRequest)
      .then((response) => {
        setLoading(false);
        setSnackbarVisible(true);
        if (response.status !== 200) {
          onFailedLogin(response.status);
          return null;
        }
        return response.json();
      })
      .then(async (data) => {
        if (data != null) {
          await onSuccessfulLogin(data);
        }
      });
  }, [setLoading, loginDetails, setErrors]);

  const boxStyle = {
    position: 'absolute',
    top: '50%',
    left: '50%',
    transform: 'translate(-50%, -50%)',
    width: 400,
    bgcolor: 'background.paper',
    boxShadow: 24,
    p: 4,
    borderRadius: 3,
  };

  const { close, visible } = props;

  return (
    <div>
      <Modal
        onClose={close}
        open={visible}
        aria-labelledby="simple-modal-title"
        aria-describedby="simple-modal-description"
      >
        <Box sx={boxStyle}>
          <Grid>
            <Grid container direction="row" justifyContent="center">
              <Grid xs item />
              <Grid align="center">
                <Avatar src={logo} sx={{ width: 72, height: 72 }} />
                <h2>Sign In</h2>
              </Grid>
              <Grid xs align="right" item>
                <CloseIcon onClick={onClose} className={styles.closeImg} />
              </Grid>
            </Grid>
            <TextField
              id="email"
              type="email"
              error={!!errors.email}
              helperText={errors.email}
              label="E-mail"
              placeholder="Enter e-mail address"
              fullWidth
              required
              onChange={handleTextChange}
            />
            <TextField
              id="password"
              error={!!errors.password}
              helperText={errors.password}
              label="Password"
              type="password"
              margin="dense"
              placeholder="Enter password"
              fullWidth
              required
              onChange={handleTextChange}
            />
            <FormControlLabel
              control={(
                <Checkbox
                  id="remember"
                  color="primary"
                  onChange={handleCheckedChange}
                />
              )}
              label="Remember me"
            />
            <LoadingButton
              loading={loading}
              type="submit"
              color="primary"
              variant="contained"
              fullWidth
              onClick={doLogin}
            >
              Sign in
            </LoadingButton>
            <Typography mt={1}>
              <Link href="#" onClick={handlepasswordForget} className={styles.pwReminderLink}>Forgot password?</Link>
            </Typography>
          </Grid>
        </Box>
      </Modal>
      <Snackbar
        open={isSnackbarVisible}
        autoHideDuration={10000}
        onClose={() => setSnackbarVisible(false)}
      >
        <Alert
          severity={errors.response ? 'error' : 'success'}
          sx={{ width: '100%' }}
          onClose={() => setSnackbarVisible(false)}
        >
          {errors.response ? errors.response : 'Login successful!'}
        </Alert>
      </Snackbar>
      <PasswordForgetForm visible={forgotPassword} onClose={() => setForgotPassword(false)} />
    </div>
  );
}

export default LoginForm;

// TODO: make the pop-up more modern and nice
