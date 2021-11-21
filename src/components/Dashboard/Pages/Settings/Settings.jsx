import React, { useCallback, useState } from 'react';
import Stack from '@mui/material/Stack';
import { Link } from 'react-router-dom';
import TextField from '@mui/material/TextField';
import Button from '@mui/material/Button';
import { createTheme, ThemeProvider, styled } from '@mui/material/styles';
import Switch from '@mui/material/Switch';
import { FormControlLabel } from '@material-ui/core';
import validator from 'validator';
import styles from './settings.module.css';
import profileDefault from '../../../../assets/profiledefault.png';

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

const IOSSwitch = styled((props) => (
  <Switch focusVisibleClassName=".Mui-focusVisible" disableRipple {...props} />
))(({ theme }) => ({
  width: 42,
  height: 26,
  padding: 0,
  '& .MuiSwitch-switchBase': {
    padding: 0,
    margin: 2,
    transitionDuration: '300ms',
    '&.Mui-checked': {
      transform: 'translateX(16px)',
      color: '#fff',
      '& + .MuiSwitch-track': {
        backgroundColor: theme.palette.mode === 'dark' ? '#1888ff' : '#1888ff',
        opacity: 1,
        border: 0,
      },
      '&.Mui-disabled + .MuiSwitch-track': {
        opacity: 0.5,
      },
    },
    '&.Mui-focusVisible .MuiSwitch-thumb': {
      color: '#1888ff',
      border: '6px solid #fff',
    },
    '&.Mui-disabled .MuiSwitch-thumb': {
      color:
                theme.palette.mode === 'light'
                  ? theme.palette.grey[100]
                  : theme.palette.grey[600],
    },
    '&.Mui-disabled + .MuiSwitch-track': {
      opacity: theme.palette.mode === 'light' ? 0.7 : 0.3,
    },
  },
  '& .MuiSwitch-thumb': {
    boxSizing: 'border-box',
    width: 22,
    height: 22,
  },
  '& .MuiSwitch-track': {
    borderRadius: 26 / 2,
    backgroundColor: theme.palette.mode === 'light' ? 'darkgray' : 'darkgray',
    opacity: 1,
    transition: theme.transitions.create(['background-color'], {
      duration: 500,
    }),
  },
}));

function PasswordSection(props) {
  if (props.changePassword) {
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
            onChange={props.onPasswordInfoChange}
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
            onChange={props.onPasswordInfoChange}
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
            onChange={props.onPasswordInfoChange}
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

  /* Input Validation */
  function isValidInput() {
    let currentErrors = {};
    if (!validator.isEmail(userInfo.email)) {
      currentErrors = { ...currentErrors, email: 'Please provide a valid email' };
    }
    if (userInfo.firstname === '') {
      currentErrors = { ...currentErrors, firstname: 'Please enter your first name!' };
    }
    if (userInfo.lastname === '') {
      currentErrors = { ...currentErrors, lastname: 'Please enter your last name!' };
    }
    setErrors(currentErrors);
    return !Object.keys(currentErrors).length;
  }

  const save = useCallback(() => {
      if (!validateInput()) {
        return;
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
              alt="profile-picture"
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
              width: '600'
                                + 'px',
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
            control={(
              <IOSSwitch
                sx={{ m: 1 }}
                check={changePassword}
                onClick={handleChangePassword}
              />
                            )}
          />
        </Stack>

        <PasswordSection
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
