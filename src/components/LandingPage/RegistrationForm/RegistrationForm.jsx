import {
  Avatar, Box, Grid, Modal, Snackbar, TextField,
} from '@mui/material';
import LoadingButton from '@mui/lab/LoadingButton';
import React, { useCallback, useState } from 'react';
import validator from 'validator';
import CloseIcon from '@mui/icons-material/Close';
import { Alert } from '@mui/lab';
import { BACKEND_ADDRESS } from '../../common/constants';
import logo from '../../../assets/logo.svg';
import styles from './registrationform.module.css';

function RegistrationForm(props) {
  const [userDetails, setUserDetails] = useState({
    email: '',
    firstname: '',
    lastname: '',
    affiliation: '',
    password: '',
    passwordAgain: '',
  });

  const [errors, setErrors] = useState({});
  const [loading, setLoading] = useState(false);
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  const handleTextChange = useCallback((e) => {
    setUserDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [setUserDetails]);

  function validateInput() {
    let currentErrors = {};
    if (!validator.isEmail(userDetails.email)) {
      currentErrors = { ...currentErrors, email: 'A valid e-mail is required!' };
    }
    if (userDetails.firstname === '') {
      currentErrors = { ...currentErrors, firstname: 'Please enter your first name!' };
    }
    if (userDetails.lastname === '') {
      currentErrors = { ...currentErrors, lastname: 'Please enter your last name!' };
    }
    if (userDetails.password === '') {
      currentErrors = { ...currentErrors, password: 'Password is required!' };
    }
    if (userDetails.passwordAgain !== userDetails.password) {
      currentErrors = { ...currentErrors, passwordAgain: 'The passwords must match!' };
    }
    setErrors(currentErrors);
    return !Object.keys(currentErrors).length;
  }

  const onClose = useCallback(() => {
    setUserDetails({
      email: '',
      firstname: '',
      lastname: '',
      affiliation: '',
      password: '',
      passwordAgain: '',
    });
    setErrors({});
    setLoading(false);
    props.onClose();
  }, [setUserDetails, setErrors, setLoading, props]);

  function onSuccessfulRegistration() {
    onClose();
    props.onSuccessfulRegistration();
  }

  function onFailedRegistration(response) {
    switch (response.status) {
      case 400:
        setErrors((prevState) => ({ ...prevState, response: 'Please check your input!' }));
        break;
      case 409:
        setErrors((prevState) => ({ ...prevState, response: 'Account already exists!' }));
        break;
      default:
        console.log(response);
        setErrors((prevState) => ({ ...prevState, response: 'Unknown error, please try again later!' }));
        break;
    }
  }

  const doRegistration = useCallback(() => {
    if (!validateInput()) {
      return;
    }
    setLoading(true);

    const requestOptions = {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        first_name: userDetails.firstname,
        last_name: userDetails.lastname,
        email: userDetails.email,
        password: userDetails.password,
        note: userDetails.affiliation,
      }),
    };
    fetch(`${BACKEND_ADDRESS}/register`, requestOptions)
      .then((response) => {
        setLoading(false);
        if (response.status === 201) {
          onSuccessfulRegistration();
        } else {
          onFailedRegistration(response);
        }
        setSnackbarVisible(true);
      });
  }, [userDetails, setErrors, setLoading, props, setSnackbarVisible]);

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
                <h2>Register new account</h2>
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
              margin="dense"
              onChange={handleTextChange}
            />
            <TextField
              id="firstname"
              type="text"
              error={!!errors.firstname}
              helperText={errors.firstname}
              label="First name"
              placeholder="Enter your first name"
              margin="dense"
              fullWidth
              required
              onChange={handleTextChange}
            />
            <TextField
              id="lastname"
              type="text"
              error={!!errors.lastname}
              helperText={errors.lastname}
              label="Last name"
              placeholder="Enter your last name"
              margin="dense"
              fullWidth
              required
              onChange={handleTextChange}
            />
            <TextField
              id="affiliation"
              error={!!errors.affiliation}
              helperText={errors.affiliation}
              label="Academic affiliation"
              type="text"
              placeholder="Enter your university, company or etc."
              margin="dense"
              fullWidth
              onChange={handleTextChange}
            />
            <TextField
              id="password"
              error={!!errors.password}
              helperText={errors.password}
              label="Password"
              type="password"
              placeholder="Enter password"
              margin="dense"
              fullWidth
              required
              onChange={handleTextChange}
            />
            <TextField
              id="passwordAgain"
              error={!!errors.passwordAgain}
              helperText={errors.passwordAgain}
              label="Password again"
              type="password"
              placeholder="Enter your password again"
              margin="dense"
              fullWidth
              required
              onChange={handleTextChange}
            />
            <Box mt={1}>
              <LoadingButton
                loading={loading}
                type="submit"
                color="primary"
                variant="contained"
                fullWidth
                onClick={doRegistration}
              >
                Sign up
              </LoadingButton>
            </Box>
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
          {errors.response ? errors.response : 'Successful registration, now you can log in!'}
        </Alert>
      </Snackbar>
    </div>
  );
}

export default RegistrationForm;
