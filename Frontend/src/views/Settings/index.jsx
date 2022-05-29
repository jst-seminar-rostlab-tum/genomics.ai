import React, { useCallback, useState } from 'react';
import Stack from '@mui/material/Stack';
import TextField from '@mui/material/TextField';
import FormControlLabel from '@mui/material/FormControlLabel';
import Switch from '@mui/material/Switch';
import { Alert, Snackbar } from '@mui/material';
import validator from 'validator';

import Button from 'components/CustomButton';
import ProfileImage from 'components/ProfileImage';

import styles from './settings.module.css';
import { useAuth } from 'shared/context/authContext';
import ProfileService from 'shared/services/Profile.service';


function PasswordSection({ onPasswordInfoChange, errors, changePassword }) {
  if (changePassword) {
    return (
      <Stack
        spacing={5}
        direction="column"
      >
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
  const [user] = useAuth();

  // Password Change
  const [changePassword, setChangePassword] = React.useState(false);
  const handleChangePassword = () => setChangePassword(!changePassword);
  // Errors
  const [errors, setErrors] = useState({});
  const [saveInProgress, setSaveInProgress] = useState(false);
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  /* Input Text Fields */
  // General Information
  const [userInfo, setUserInfo] = useState({
    firstName: user.firstName,
    lastName: user.lastName,
    emailAddress: user.email,
    academicAffiliation: user.note,
  });

  const onUserInfoChange = useCallback((e) => {
    setUserInfo((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [userInfo]);

  // Password
  const [passwordInfo, setPasswordInfo] = useState({
    newPassword: '',
    newPasswordRepeated: '',
  });

  const onPasswordInfoChange = (e) => {
    setPasswordInfo((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  };

  /* Input Validation */
  function isValidInput() {
    let currentErrors = {};
    if (!validator.isEmail(userInfo.emailAddress)) {
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

  function isValidPasswordInfo() {
    let newErrors = {};
    // password correctly repeated?
    if (passwordInfo.newPassword !== passwordInfo.newPasswordRepeated) {
      newErrors = { ...newErrors, newPasswordRepeated: 'The passwords must match!' };
    }
    setErrors((currentErrors) => ({ ...currentErrors, ...newErrors }));
    console.log(`${passwordInfo.newPassword} ${passwordInfo.newPasswordRepeated}`);
    return !Object.keys(newErrors).length;
  }

  const saveUserData = () => {
    if (isValidInput() && isValidPasswordInfo()) {
      setSaveInProgress(true);
      ProfileService.updateProfile(userInfo, changePassword ? passwordInfo.newPassword : null)
        .catch((err) => {
          setErrors((prevState) => ({ ...prevState, response: "Couldn't save changes :/" }));
          console.log(err);
        }).finally(() => {
          setSaveInProgress(false);
          setSnackbarVisible(true);
        });
    }
  };
  return (
    <>
      <Stack
        direction="column"
        sx={{
          paddingTop: '100px',
          paddingLeft: '130px',
          display: 'flex',
          width: '100%',
          justifyContent: 'center',
        }}
      >
        <div className={styles.title}>
          <h1>Settings</h1>
        </div>

        <Stack
          className={styles.background}
          spacing={5}
          direction="column"
        >
          {/* Profile Image --------------------------------------------------------------*/}
          <div className={styles.profilePicture}>
            <ProfileImage sizePixels={180} editable />
          </div>

          {/* Input Fields --------------------------------------------------------------*/}
          <Stack
            spacing={2}
            style={{ display: 'flex' }}
            direction="row"
          >
            <div className={styles.inputComponent}>
              <div className={styles.inputText}> First Name </div>
              <TextField
                id="firstName"
                value={userInfo.firstName}
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
                value={userInfo.lastName}
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
              disabled
              value={userInfo.emailAddress}
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
            <div className={styles.inputText}>Academic Affiliation</div>
            <TextField
              id="academicAffiliation"
              value={userInfo.academicAffiliation}
              label="Academic Affiliation (e.g. University)"
              type="text"
              style={{ width: '600px', background: 'white' }}
              onChange={onUserInfoChange}
            />
          </div>

          <div className={styles.inputComponent}>
            <div className={styles.inputText}>Password</div>
            <FormControlLabel
              label="Change Password"
              sx={{ m: 1 }}
              control={(
                <Switch
                  checked={changePassword}
                  onClick={handleChangePassword}
                />
              )}
            />
          </div>

          <PasswordSection
            errors={errors}
            changePassword={changePassword}
            passwordInfo={passwordInfo}
            onPasswordInfoChange={onPasswordInfoChange}
          />

          <div style={{
            display: 'flex',
            justifyContent: 'flex-start',
          }}
          >
            <Button disabled={saveInProgress} onClick={saveUserData}>
              {saveInProgress ? 'Saving...' : 'Save'}
            </Button>
          </div>
        </Stack>

      </Stack>

      {/* Snackbar --------------------------------------------------------------*/}
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
          {errors.response ? errors.response : 'Changes saved.'}
        </Alert>
      </Snackbar>
    </>
  );
}

export default Settings;
