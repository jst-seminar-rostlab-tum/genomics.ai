/* eslint-disable */
import React, { useCallback, useState } from "react";
import {
  Modal,
  Box,
  Grid,
  TextField,
  Button,
  Avatar,
  Snackbar,
  Alert,
  Typography,
} from "@mui/material";
import CloseIcon from "@mui/icons-material/Close";
import { useHistory } from "react-router-dom";
import validator from "validator";
import { BACKEND_ADDRESS } from "../../shared/utils/common/constants";
import styles from "./passwordforgetform.module.css";
import logo from "../../assets/logo.svg";
import Input from "../Input/Input";
import CustomButton from "components/CustomButton";

function PasswordForgetForm(props) {
  const { close, visible, switchForm } = props;
  const [errors, setErrors] = useState({});
  // eslint-disable-next-line no-unused-vars
  const [email, setEmail] = useState();
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  const history = useHistory();

  const boxStyle = {
    position: "absolute",
    top: "50%",
    left: "50%",
    transform: "translate(-50%, -50%)",
    width: 400,
    bgcolor: "background.paper",
    boxShadow: 24,
    p: 4,
    borderRadius: 3,
  };

  const onClose = useCallback(() => {
    setEmail("");
    setErrors({});
    props.onClose();
  }, [props]);

  function onSendClick() {
    if (!validateInput()) {
      return;
    }

    const requestOptions = {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        email,
      }),
    };
    fetch(`${BACKEND_ADDRESS}/password_reset`, requestOptions).then(
      (response) => {
        if (response.status == 404) {
          setErrors({ response: "User not found" });
        }
        if (response.status == 400) {
          setErrors({ response: response.statusText });
        } else {
        }
        setSnackbarVisible(true);
      }
    );
    
    //props.onClick()
    //switchForm(true);

    // history.push('/password_reset');
  }

  function onSnackbarClose() {
    setSnackbarVisible(false);
    setErrors({});
  }

  function validateInput() {
    if (email == null) {
      setErrors({ ...errors, email: "Empty input!" });
      return false;
    }

    if (!validator.isEmail(email)) {
      setErrors({ ...errors, email: "A valid e-mail is required!" });
      return false;
    }
    return true;
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
                <h2>Forgot Password</h2>
              </Grid>
              <Grid xs align="right" item>
                <CloseIcon onClick={onClose} className={styles.closeImg} />
              </Grid>
            </Grid>
            <Input
              id="email"
              type="email"
              error={!!errors.email}
              helperText={errors.email}
              label="E-mail"
              placeholder="Enter e-mail address"
              isRequired={true}
              onChangeEvent={(e) => setEmail(e.target.value)}
            />
            <Box
              sx={{ width: "100%", display: "flex", justifyContent: "center"}}
            >
              <CustomButton type="primary" onClick={onSendClick}>
                <Typography>Send</Typography>
              </CustomButton>
            </Box>
          </Grid>
        </Box>
      </Modal>

      <Snackbar
        open={isSnackbarVisible}
        autoHideDuration={5000}
        onClose={onSnackbarClose}
      >
        <Alert
          severity={errors.response ? "error" : "success"}
          sx={{ width: "100%" }}
          onClose={onSnackbarClose}
        >
          {errors.response
            ? errors.response
            : "Password reset email has been sent! Please check your e-mail!"}
        </Alert>
      </Snackbar>
    </div>
  );
}

export default PasswordForgetForm;