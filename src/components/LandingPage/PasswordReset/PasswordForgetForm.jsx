import React, { useCallback, useState } from 'react';
import {
  Modal, Box, Grid, TextField, Button, Avatar, Snackbar,
} from '@mui/material';
import CloseIcon from '@mui/icons-material/Close';
import { useHistory } from 'react-router-dom';
import { Alert } from '@mui/lab';
import validator from 'validator';
import styles from './passwordreset.module.css';
import logo from '../../../assets/logo.svg';

function PasswordForgetForm(props) {
  const { close, visible } = props;
  const [errors, setErrors] = useState({});
  // eslint-disable-next-line no-unused-vars
  const [email, setEmail] = useState();
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  const history = useHistory();

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

  const onClose = useCallback(() => {
    setEmail('');
    setErrors({});
    props.onClose();
  }, [setEmail, setErrors, props]);

  function onSendClick () {
    setSnackbarVisible(true);
    if(!validateInput()){
      return;
    }
    
    history.push('/passwordreset');
  }

  function validateInput() {
    if (!validator.isEmail(email)) {
      setErrors({ ...errors,email : 'A valid e-mail is required!' });
      return false;
    }
    else return true;
  }


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
                <h2>Forget Password</h2>
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
              margin="dense"
              required
              fullWidth
              onChange={() => setEmail()}
            />
            <Button
              type="submit"
              variant="outlined"
              fullWidth
              sx={{
                pt: 1,
              }}
              onClick={onSendClick}
            >
              send confirm
            </Button>

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
          {errors.response ? errors.response : 'Password reset email has been sent'}
        </Alert>
      </Snackbar>        

    </div>
  );
}

export default PasswordForgetForm;
