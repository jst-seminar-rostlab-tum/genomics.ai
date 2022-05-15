import {
  Avatar, Typography, Box, Grid, Snackbar, TextField,
} from '@mui/material';
import LoadingButton from '@mui/lab/LoadingButton';
import React, { useCallback, useState } from 'react';
import validator from 'validator';
import CloseIcon from '@mui/icons-material/Close';
import { Alert } from '@mui/lab';
import { BACKEND_ADDRESS } from 'shared/utils/common/constants';
import logo from 'assets/logo.svg';
import styles from './registrationform.module.css';
import { useAuth } from 'shared/context/authContext';
import Input from 'components/Input/Input';
import { Modal, ModalTitle } from 'components/Modal';
import CustomButton from 'components/CustomButton';

function RegistrationForm(props) {
  const [, setUser] = useAuth();
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
  }, [setUserDetails, setErrors, setLoading, props, setUser]);

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
        if (response.status === 200 || response.status === 201) {
          onSuccessfulRegistration();
        } else {
          onFailedRegistration(response);
        }
        setSnackbarVisible(true);
      });
  }, [userDetails, setErrors, setLoading, props, setSnackbarVisible, setUser]);

  const { close, visible } = props;

  return (
    <div>
      <Modal
        setOpen={(o) => !o && onClose()}
        isOpen={visible}
      >
        <ModalTitle>Register new user</ModalTitle>
        <Box sx={{ width: 340 }}>
          <Grid sx={{ width: 320 }}>
            <Input
              id="email"
              type="email"
              error={!!errors.email}
              helperText={errors.email}
              label="E-mail"
              placeholder="Enter e-mail address"
              isRequired
              onChangeEvent={handleTextChange}
            />
            <Input
              id="firstname"
              type="text"
              error={!!errors.firstname}
              helperText={errors.firstname}
              label="First name"
              placeholder="Enter your first name"
              isRequired
              onChangeEvent={handleTextChange}
            />
            <Input
              id="lastname"
              type="text"
              error={!!errors.lastname}
              helperText={errors.lastname}
              label="Last name"
              placeholder="Enter your last name"
              isRequired
              onChangeEvent={handleTextChange}
            />
            <Input
              id="affiliation"
              error={!!errors.affiliation}
              helperText={errors.affiliation}
              label="Academic affiliation"
              type="text"
              placeholder="Enter your university, company or etc."
              onChangeEvent={handleTextChange}
            />
            <Input
              id="password"
              error={!!errors.password}
              helperText={errors.password}
              label="Password"
              type="password"
              placeholder="Enter password"
              isRequired
              onChangeEvent={handleTextChange}
            />
            <Input
              id="passwordAgain"
              error={!!errors.passwordAgain}
              helperText={errors.passwordAgain}
              label="Password again"
              type="password"
              placeholder="Enter your password again"
              isRequired
              onChangeEvent={handleTextChange}
            />
            <Box mt={1}>
              <CustomButton type="contained" sx={{ mr: 1 }} onClick={doRegistration}>
                <Typography>Sign up</Typography>
              </CustomButton>
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
          {errors.response ? errors.response : 'Successful registration, check your e-mails for verification!'}
        </Alert>
      </Snackbar>
    </div>
  );
}

export default RegistrationForm;

// TODO: make the popup more modern and nice
