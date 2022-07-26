/* eslint-disable */
import {
  Alert,
  Box,
  Checkbox,
  FormControlLabel,
  Grid,
  Snackbar,
  Typography,
  Link,
} from "@mui/material";
import React, { useCallback, useState, useEffect, useContext } from "react";
import validator from "validator";
import CloseIcon from "@mui/icons-material/Close";
import { useHistory } from "react-router-dom";
import logo from "assets/logo.svg";
import styles from "./loginform.module.css";
import { BACKEND_ADDRESS } from "shared/utils/common/constants";
import PasswordForgetForm from "../PasswordForgetForm";
import { useAuth } from "shared/context/authContext";
import { Modal, ModalTitle } from "components/Modal";
import Input from "components/Input/Input";
import CustomButton from "components/CustomButton";
import { colors } from 'shared/theme/colors';
import { LoginContext } from "shared/context/loginContext";

function LoginForm(props) {
  const [, setUser] = useAuth();

  const [loginDetails, setLoginDetails] = useState({
    email: "",
    password: "",
    remember: false,
  });

  const [errors, setErrors] = useState({});
  const [loading, setLoading] = useState(false);
  const [isSnackbarVisible, setSnackbarVisible] = useState(false);

  const { loginClose: onClose, loginVisible: visible, switchRegister: switchForm, switchForget: onForgetPassword } = useContext(LoginContext);

  const history = useHistory();

  useEffect(() => {
    if(localStorage.getItem("email")&&localStorage.getItem("email").length !==0) {
      setLoginDetails({email:localStorage.getItem("email")})
    }
  },[])
  
  const clearForm = () => {
    setLoginDetails({
      email: "",
      password: "",
      remember: false,
    });
    setErrors({});
    setLoading(false);
    onClose();
  }

  const handleTextChange = useCallback(
    (e) => {
      setLoginDetails((prevState) => ({
        ...prevState,
        [e.target.id]: e.target.value,
      }));
    },
    [setLoginDetails]
  );

  const handleCheckedChange = useCallback(
    (e) => {
      setLoginDetails((prevState) => ({
        ...prevState,
        [e.target.id]: e.target.checked,
      }));
    },
    [setLoginDetails]
  );

  function validateInput() {
    let currentErrors = {};
    if (!validator.isEmail(loginDetails.email)) {
      currentErrors = {
        ...currentErrors,
        email: "A valid e-mail is required!",
      };
    }
    if (loginDetails.password === "") {
      currentErrors = { ...currentErrors, password: "Password is required!" };
    }
    setErrors(currentErrors);
    return !Object.keys(currentErrors).length;
  }

  async function onSuccessfulLogin(data) {
    localStorage.setItem("jwt", data.jwt);
    localStorage.setItem("user", JSON.stringify(data.user));
    if(loginDetails.remember) {
      localStorage.setItem("email",loginDetails.email);
    }
    onClose();
    history.go(0);
    await setUser(data.user);
  }

  function onFailedLogin(code) {
    switch (code) {
      case 401:
        setErrors((prevState) => ({
          ...prevState,
          response: "Wrong credentials!",
        }));
        break;
      case 404:
        setErrors((prevState) => ({
          ...prevState,
          response: "User not found!",
        }));
        break;
      default:
        setErrors((prevState) => ({
          ...prevState,
          response: "Unknown error, please try again later!",
        }));
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
      method: "POST",
      headers: {
        "Content-Type": "application/json",
        Authorization: "Bearer",
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

    history.push("/");
  }, [setLoading, loginDetails, setErrors]);


  
  return (
    <div>
      <Modal setOpen={(o) => !o && onClose()} isOpen={visible}>
        <ModalTitle>Log in</ModalTitle>
        <Box>
          <Grid sx={{ pr: "1em" }}>
            <Input
              id="email"
              type="email"
              defaultValue={localStorage.getItem("email")}
              
              error={!!errors.email}
              helperText={errors.email}
              label="E-mail"
              placeholder="Enter e-mail address"
              isRequired={true}
              onChangeEvent={handleTextChange}
            />
            <Input
              id="password"
              error={!!errors.password}
              helperText={errors.password}
              label="Password"
              type="password"
              placeholder="Enter password"
              isRequired={true}
              onChangeEvent={handleTextChange}
            />
            <Box sx={{ paddingLeft: 1, display: "flex", flexDirection: "column", gap: "4   px" }}>
              <FormControlLabel
                control={
                  <Checkbox
                    id="remember"
                    sx={{ color: colors.primary[400], '& .MuiSvgIcon-root': { color: colors.primary[400] } }}
                    onChange={handleCheckedChange}
                    disableRipple
                    disableFocusRipple
                  />
                }
                label="Remember me"
              />
              <CustomButton type="primary" onClick={doLogin}>
                <Typography>Sign in</Typography>
              </CustomButton>
              <Typography mt={1} textAlign="center">
                <Link
                  href="#"
                  onClick={() => { onForgetPassword(true); clearForm(); }}
                >
                  Forgot password?
                </Link>
              </Typography>
              <Typography mt={1} textAlign="center">
                Don't have an account yet? Sign up {" "}
                <Link
                  href="#"
                  onClick={() => { switchForm(true); clearForm(); }} 
                > 
                  here
                </Link>
                {"."}
              </Typography>
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
          severity={errors.response ? "error" : "success"}
          sx={{ width: "100%" }}
          onClose={() => setSnackbarVisible(false)}
        >
          {errors.response ? errors.response : "Login successful!"}
        </Alert>
      </Snackbar>
    </div>
  );
}

export default LoginForm;