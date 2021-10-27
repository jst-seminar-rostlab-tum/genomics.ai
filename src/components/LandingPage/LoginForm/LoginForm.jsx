import {Avatar, Box, Checkbox, FormControlLabel, Grid, Modal, TextField, Typography} from "@mui/material";
import LoadingButton from '@mui/lab/LoadingButton';
import logo from "../../../assets/logo.svg"
import {useState} from "react";
import validator from 'validator'
import styles from './loginform.module.css'
import CloseIcon from '@mui/icons-material/Close';
import { BACKEND_ADDRESS } from "../../common/constants";

function LoginForm(props) {
    const [loginDetails, setLoginDetails] = useState({
        email: "",
        password: "",
        remember: false,
    })

    const [errors, setErrors] = useState({})
    const [loading, setLoading] = useState(false)
    const [jwt , setJwt] = useState("")
    

    function handleTextChange(e) {
        setLoginDetails(prevState => ({...prevState, [e.target.id]: e.target.value}));
    }

    function handleCheckedChange(e) {
        setLoginDetails(prevState => ({...prevState, [e.target.id]: e.target.checked}));
    }

    function validateInput() {
        let errors = {};
        if (!validator.isEmail(loginDetails.email)) {
            errors = {...errors, email: "A valid e-mail is required!"};
        }
        if (loginDetails.password === "") {
            errors = {...errors, password: "Password is required!"};
        }
        setErrors(errors);
        return !Object.keys(errors).length;
    }

    function doLogin() {
        if (!validateInput()) {
            return;
        }
        setLoading(true);
        //server communication 

        const loginRequest = {
            method: 'POST',
            headers: {
            'Content-Type': 'application/json',
            'Authorization': 'Bearer',
            },
            body: JSON.stringify({
                email: loginDetails.email,
                password: loginDetails.password,
            })
        };
        fetch(BACKEND_ADDRESS + "/auth", loginRequest)
            .then(response => {
                setLoading(false);
                if (response.status === 200) {
                onSuccessfulLogin();
                } else {
                    onFailedLogin(response.status);
                    return;
                }
                return response.json();
            })
            .then(data => { 
                if(data != null){
                    setJwt(data.jwt);
                    localStorage.jwt = data.jwt;
                    console.log(data.msg);
                }
                
            })
    }


    function onSuccessfulLogin(){
        onClose();
        /////TODO:switch to dashboard////
    }

    function onFailedLogin(code){
        switch (code) {
            case 401:
              setErrors(prevState => ({...prevState, response: "wrong credentials"}));
              break;
            case 404:
              setErrors(prevState => ({...prevState, response: "user not found"}));
              break;
            default:
              setErrors(prevState => ({...prevState, response: "Unknown error, please try again later!"}));
              break;
          }
    }

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

    function onClose() {
        setLoginDetails({
            email: "",
            password: "",
            remember: false,
        });
        setErrors({});
        setLoading(false);
        props.onClose();
    }

    return (
        <div>
            <Modal
                onClose={props.onClose}
                open={props.visible}
                aria-labelledby="simple-modal-title"
                aria-describedby="simple-modal-description">
                <Box sx={boxStyle}>
                    <Grid>
                        <Grid container direction={'row'} justifyContent={'center'}>
                            <Grid xs item/>
                            <Grid align={'center'}>
                                <Avatar src={logo} sx={{width: 72, height: 72}}/>
                                <h2>Sign In</h2>
                            </Grid>
                            <Grid xs align={'right'} item>
                                <CloseIcon onClick={onClose} className={styles.closeImg}/>
                            </Grid>
                        </Grid>
                        <TextField id='email'
                                   type='email'
                                   error={!!errors['email']}
                                   helperText={errors['email']}
                                   label='E-mail'
                                   placeholder='Enter e-mail address'
                                   fullWidth required
                                   onChange={handleTextChange}/>
                        <TextField id='password'
                                   error={!!errors['password']}
                                   helperText={errors['password']}
                                   label='Password'
                                   type='password'
                                   margin='dense'
                                   placeholder='Enter password'
                                   fullWidth required
                                   onChange={handleTextChange}/>
                        <FormControlLabel
                            control={
                                <Checkbox
                                    id="remember"
                                    color="primary"
                                    onChange={handleCheckedChange}
                                />
                            }
                            label="Remember me"
                        />
                        <LoadingButton loading={loading}
                                       type='submit'
                                       color='primary'
                                       variant="contained"
                                       fullWidth
                                       onClick={doLogin}>Sign in
                        </LoadingButton>
                        <Typography mt={1}>
                            <a href="/" className={styles.pwReminderLink}>
                                Forgot password?
                            </a>
                        </Typography>
                    </Grid>
                </Box>
            </Modal>
            



        </div>
    );
}

export default LoginForm;