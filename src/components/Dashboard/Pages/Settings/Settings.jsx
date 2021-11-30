import React, { useCallback, useState } from 'react';
import Stack from '@mui/material/Stack';
import { Link } from 'react-router-dom';
import TextField from '@mui/material/TextField';
import Button from '@mui/material/Button';
import FormControlLabel from '@mui/material/FormControlLabel';
import Switch from '@mui/material/Switch';
import { createTheme, ThemeProvider } from '@mui/material/styles';
import validator from 'validator';
import styles from './settings.module.css';
import profileDefault from '../../../../assets/user.png';

const myTheme = createTheme({
  palette: {
    primary: {
      main: '#1888ff',
    },
    secondary: {
      main: '#6B7379',
    },
  },
});

function PasswordSection({ onPasswordInfoChange, errors, changePassword }) {
  if (changePassword) {
    return (
      <Stack
        spacing={5}
        direction="column"
      >

        <div className={styles.inputComponent}>
          <div className={styles.inputText}> Current Password </div>
          <TextField
            id="currentPassword"
            required
            label="Current Password"
            type="password"
            style={{ width: '600px', background: 'white' }}
            onChange={onPasswordInfoChange}
            error={!!errors.currentPassword}
            helperText={errors.currentPassword}
          />
        </div>

        <div className={styles.inputComponent}>
          <div className={styles.inputText}> New password </div>
          <TextField
            id="newPassword"
            required
            label="New Password"
            type="password"
            style={{ width: '600px', background: 'white' }}
            onChange={onPasswordInfoChange}
            error={!!errors.newPassword}
            helperText={errors.newPassword}
          />
        </div>

        <div className={styles.inputComponent}>
          <div className={styles.inputText}> Repeat new password </div>
          <TextField
            id="newPasswordRepeated"
            required
            label="Repeat old password"
            type="password"
            style={{ width: '600px', background: 'white' }}
            onChange={onPasswordInfoChange}
            error={!!errors.newPasswordRepeated}
            helperText={errors.newPasswordRepeated}
          />
        </div>
      </Stack>
    );
  }
  return null;
}

function Settings() {
  /* Booleans */
  // Password Change
  const [changePassword, setChangePassword] = React.useState(false);
  const handleChangePassword = () => setChangePassword(!changePassword);
  // Errors
  const [errors, setErrors] = useState({});

  /* Input Text Fields */
  // General Information
  const [userInfo, setUserInfo] = useState({
    firstName: '',
    lastName: '',
    emailAddress: '',
    academicAffiliation: '',
    aboutMe: '',
  });

  const onUserInfoChange = useCallback((e) => {
    setUserInfo((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [userInfo]);

  // Password
  const [passwordInfo, setPasswordInfo] = useState({
    currentPassword: '',
    newPassword: '',
    newPasswordRepeated: '',
  });

  const onPasswordInfoChange = useCallback((e) => {
    setPasswordInfo((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [passwordInfo]);

  function passwordIsSafe(minLength) {
    // minimum length
    if (userInfo.newPassword.length < minLength) {
      return 'Password must be at least 8 characters long!';
    }
    // password should contain these
    let lowerCase = false;
    let upperCase = false;
    let numbers = false;

    userInfo.newPassword.split('').forEach((c) => {
      if (c === c.toUpperCase()) {
        upperCase = true;
      } else if (c === c.toLowerCase()) {
        lowerCase = true;
      } else if (!Number.isNaN(c)) {
        numbers = true;
      }
    });

    if (!lowerCase || !upperCase || !numbers) {
      return 'Use numbers, upper & lower case letters!';
    }
    return '';
  }

  /* Input Validation */
  function isValidInput() {
    let currentErrors = {};
    if (!validator.isEmail(userInfo.email)) {
      currentErrors = { ...currentErrors, emailAddress: 'Please provide a valid email' };
    }
    if (userInfo.firstName === '') {
      currentErrors = { ...currentErrors, firstname: 'Please enter your first name!' };
    }
    if (userInfo.lastName === '') {
      currentErrors = { ...currentErrors, lastname: 'Please enter your last name!' };
    }
    setErrors(currentErrors);
    return !Object.keys(currentErrors).length;
  }

  function currentPasswordIsCorrect() {
    // validate password
    return true;
  }

  function isValidPasswordInfo() {
    let newErrors = {};
    // current password correct?
    if (!currentPasswordIsCorrect()) {
      newErrors = { ...newErrors, currentPassword: 'Please enter your last name!' };
    }
    // password correctly repeated?
    if (userInfo.newPassword === userInfo.newPasswordRepeated) {
      newErrors = { ...newErrors, newPasswordRepeated: 'The passwords must match!' };
    }
    // new password save enough?
    const safetyErr = passwordIsSafe(8);
    if (safetyErr !== '') {
      newErrors = { ...newErrors, newPassword: safetyErr };
    }
    setErrors((currentErrors) => [...currentErrors, ...newErrors]);
    return !Object.keys(newErrors).length;
  }

  const saveUserData = useCallback(() => {
    if (isValidInput() && isValidPasswordInfo()) {
    }
  }, []);

  return (
    <>
      <Stack
        className={styles.background}
        style={{ display: 'flex' }}
        spacing={5}
        direction="column"
      >

        <div style={{ paddingBlock: '25px', paddingLeft: '80px', marginBottom: '50px' }}>
          <Link
            className={styles.profilePicture}
            to="#"
          >
            <img
              alt="profile"
              src={profileDefault}
              style={{ height: '150px' }}
            />
          </Link>
        </div>

        <Stack
          spacing={2}
          style={{ display: 'flex' }}
          direction="row"
        >
          <div className={styles.inputComponent}>
            <div className={styles.inputText}> First Name </div>
            <TextField
              id="firstName"
              required
              label="First Name"
              type="text"
              style={{ width: '290px', background: 'white' }}
              onChange={onUserInfoChange}
              error={!!errors.firstName}
              helperText={errors.firstName}
            />
          </div>

          <div className={styles.inputComponent}>
            <div className={styles.inputText}> Last Name</div>
            <TextField
              id="lastName"
              required
              label="Last Name"
              type="text"
              style={{ width: '290px', background: 'white' }}
              onChange={onUserInfoChange}
              error={!!errors.lastName}
              helperText={errors.lastName}
            />
          </div>
        </Stack>

        <div className={styles.inputComponent}>
          <div className={styles.inputText}> Email Address </div>
          <TextField
            id="emailAddress"
            required
            label="Email"
            type="text"
            style={{ width: '600px', background: 'white' }}
            onChange={onUserInfoChange}
            error={!!errors.emailAddress}
            helperText={errors.emailAddress}
          />
        </div>

        <div className={styles.inputComponent}>
          <div className={styles.inputText}> Academic Affiliation </div>
          <TextField
            id="academicAffiliation"
            label="Academic Affiliation (e.g. University)"
            type="text"
            style={{ width: '600px', background: 'white' }}
            onChange={onUserInfoChange}
          />
        </div>

        <div className={styles.inputComponent}>
          <div className={styles.inputText}> About Me </div>
          <TextField
            id="aboutMe"
            label="Tell us something about you"
            multiline
            rows="5"
            type="text"
            style={{
              width: '600px',
              background: 'white',
            }}
            onChange={onUserInfoChange}
          />
        </div>

        <div
          className={styles.divider}
          style={{ marginTop: '40px' }}
        />

        <Stack
          spacing={5}
          direction="row"
        >
          <div className={styles.headline}> Password </div>
          <FormControlLabel
            label="Change Password"
            sx={{ m: 1 }}
            control={(
              <Switch
                check={changePassword}
                onClick={handleChangePassword}
              />
            )}
          />
        </Stack>

        <PasswordSection
          errors={errors}
          changePassword={changePassword}
          passwordInfo={passwordInfo}
          onPasswordInfoChange={onPasswordInfoChange}
        />

        <div
          className={styles.divider}
          style={{ marginTop: '40px' }}
        />

        <ThemeProvider theme={myTheme}>
          <Stack
            spacing={2}
            direction="row"
            style={{ height: '55px', width: '400px' }}
          >
            <Button
              variant="contained"
              color="primary"
              style={{ width: '185px', borderRadius: '10px', fontWeight: 'bold' }}
              onClick={saveUserData}
            >
              Save
            </Button>
          </Stack>
        </ThemeProvider>
      </Stack>
    </>
  );
}

export default Settings;
