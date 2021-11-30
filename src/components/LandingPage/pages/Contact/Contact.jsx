import React from 'react';
import NavBar from '../../../NavBar/NavBar';
import styles from './contact.module.css';
import Footer from '../../Footer/Footer';
import { TextField, Typography , Grid, Button , Box} from '@mui/material';
import { useCallback, useState } from 'react';

function Contact() {
  const [contactDetails, setContactDetails] = useState({
    email: '',
    firstname: '',
    lastname: '',
    message:''
  });

  const handleTextChange = useCallback((e) => {
    setContactDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [setContactDetails]);


  const resetForm = useCallback(() => {
    setContactDetails({
      email: '',
      firstname: '',
      lastname: '',
      message: ''
    });
  }, [setContactDetails]);




  const formSubmit = async (e) => {
    e.preventDefault();
    
    let data = {
      firstname: contactDetails.firstname,
      lastname: contactDetails.lastname,
      email: contactDetails.email,
      message: contactDetails.message
    };

    try {
      // send POST to the server
      // await axios.post("BACKEND_URL", data);
      setContactDetails({
        email: '',
        firstname: '',
        lastname: '',
        message: '...sending'
      })
      
    } catch (error) {
      console.log(error);
    }
  };

  function onSuccessfulSubmit(){
    resetForm();
  }
  function onFailSubmit(){

  }

  return (
      <div className={styles.headerContainer}>
          <NavBar />
          <Typography sx={{ fontWeight: '400', fontSize: '24px' }}>Contact Us</Typography>
          <text>Write a few lines about each one and contact us about any further
                    collaboration. We will responde get back to you in a couple of
                    hours.
          </text>
            <Box  
                  component="span" 
                  margin ="auto"
                  className={styles.formContainer}
                  sx={{
                    width: "1000px",
                    height: "500px",
                    maxWidth : '100%',
                    justifyContent : "center",
                    alignItems : 'center',
                  }}  
            >
              <Grid container 
                    spacing={1} 
                    direction ="column"
                    justifyContent="center"
                    sx = {{ align: 'center',
                          display: "flex",
                    }}
              >
              
                <Grid item >
                      <TextField
                      id="email"
                      label="Email"
                      placeholder="Enter your email address"
                      variant="filled"
                      value={contactDetails.email}
                      onChange={handleTextChange}
                      required
                      fullWidth
                      type="email"
                      />
                </Grid>
                
                    <Grid item >
                      <Grid 
                          container 
                          spacing={1}
                          direction="row"
                          sx = {{ align: 'center',
                              display: "flex"
                          }}
                      >
                        <Grid 
                          item  
                          xs={6} 
                        >
                          <TextField
                          id="firstname"
                          label="Firstname"
                          placeholder="Enter your firstname"
                          variant="filled"
                          rowsMax={1}
                          value={contactDetails.firstname}
                          onChange={handleTextChange}
                          required
                          fullWidth
                          type="text"
                          />
                        </Grid>
                        <Grid 
                          item 
                          xs={6}
                          >
                          <TextField
                          id="lastname"
                          label="Lastname"
                          placeholder="Enter your lastname"
                          variant="filled"
                          rowsMax={1}
                          value={contactDetails.lastname}
                          onChange={handleTextChange}
                          fullWidth
                          required
                          type="text"
                          />
                        </Grid>

                      </Grid>
                    </Grid>

                <Grid item  >
                  <TextField
                  id="message"
                  label="Message"
                  placeholder="Enter your message"
                  variant="filled"
                  multiline 
                  rows = {8}
                  rowsMax={20}
                  value={contactDetails.message}
                  onChange={handleTextChange}
                  fullWidth
                  required
                  type="text"
                  />
                </Grid>
                <Grid item >
                    <Button type="submit" variant="outlined" color="primary" fullwidth onClick={formSubmit}> submit</Button>
                </Grid>

              </Grid>
            </Box>
            <br/>
            <br/>
            <br/>
            <br/>
            <br/>
            <br/>
            <br/>
            <br/>
            <br/>
            <br/>
            
          <Footer/>
      </div>
  )
}

export default Contact;