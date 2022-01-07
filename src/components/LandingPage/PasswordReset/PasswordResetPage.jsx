/* eslint-disable no-unused-vars */
import React, { useCallback, useState } from 'react';
import {
  Box, Grid, TextField, Button, Avatar,
} from '@mui/material';
import { useHistory } from 'react-router-dom';
import NavBar from '../../NavBar/NavBar';
import Footer from '../Footer/Footer';
import styles from './passwordreset.module.css';
import logo from '../../../assets/logo.svg';

function PasswordResetPage(props) {
  const [errors, setErrors] = useState({});
  const [newPass, setNewPass] = useState({});
  const [confirmPass, setConfirmPass] = useState({});

  const boxStyle = {
    position: 'relative',
    align: 'center',
    width: 1000,
    bgcolor: 'background.paper',
    boxShadow: 3,
    p: 12,
    borderRadius: 3,
  };

  const history = useHistory();

  const { close, visible } = props;
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
            <TextField
              id="new-password"
              type="password"
              error={!!errors.email}
              helperText={errors.email}
              label="New Password"
              placeholder="Enter New Password"
              margin="dense"
              required
              fullWidth
              onChange={setNewPass}
            />
            <TextField
              id="confirm-password"
              type="password"
              error={!!errors.email}
              helperText={errors.email}
              label="Confirm Password"
              placeholder="Reenter Password"
              margin="dense"
              required
              fullWidth
              onChange={setConfirmPass}
            />
            <Button
              type="submit"
              variant="outlined"
              fullWidth
              sx={{
                pt: 1,
              }}
              onClick={() => { history.push('/passwordreset'); }}
              size="large"
            >
              send confirm
            </Button>

          </Grid>
        </Box>
      </Box>

      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <br />
      <Footer />
    </div>
  );
}

export default PasswordResetPage;
