import React, { useCallback, useState } from 'react';
import {
  TextField, Typography, Grid, Button, Box, createTheme, ThemeProvider,
} from '@mui/material';
import Stack from '@mui/material/Stack';
import NavBar from 'components/NavBar/old';
import styles from './contact.module.css';
import Footer from 'components/Footer/old';

const theme = createTheme({
  palette: {
    primary: {
      main: '#0075FF',
    },
    secondary: {
      main: '#FFFFFF',
    },
  },
});

const Contact = () => {
  const [contactDetails, setContactDetails] = useState({
    email: '',
    firstname: '',
    lastname: '',
    message: '',
  });

  const handleTextChange = useCallback((e) => {
    setContactDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [setContactDetails]);

  const formSubmit = async (e) => {
    e.preventDefault();

    const data = {
      firstname: contactDetails.firstname,
      lastname: contactDetails.lastname,
      email: contactDetails.email,
      message: contactDetails.message,
    };
    console.log(data);

    try {
      // send POST to the server
      // await axios.post("BACKEND_URL", data);
      setContactDetails({
        email: '',
        firstname: '',
        lastname: '',
        message: '...sending',
      });
    } catch (error) {
      console.log(error);
    }
  };

  return (

    <ThemeProvider theme={theme}>
      <div>
        <NavBar />

        <div className={styles.headerContainer}>
          <Stack
            spacing="30px"
          >
            <Typography sx={{ fontWeight: 'bold', fontSize: '30px' }}>Contact Us</Typography>
            <div className={styles.textContainer}>
              <Typography
                sx={{ fontSize: '25px', paddingInline: '350px', maxWidth: '1700px' }}
                align="center"
              >
                Please message us in case you have any questions,
                feedback or collaboration-related inquiries concerning Genomics.ai.
              </Typography>
            </div>
          </Stack>
          <Box
            component="span"
            margin="auto"
            className={styles.formContainer}
            sx={{
              width: '1000px',
              height: '500px',
              maxWidth: '100%',
              justifyContent: 'center',
              alignItems: 'center',
            }}
          >
            <Grid
              container
              spacing={1}
              direction="column"
              justifyContent="center"
              sx={{
                align: 'center',
                display: 'flex',
              }}
            >

              <Grid item>
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

              <Grid item>
                <Grid
                  container
                  spacing={1}
                  direction="row"
                  sx={{
                    align: 'center',
                    display: 'flex',
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

              <Grid item>
                <TextField
                  id="message"
                  label="Message"
                  placeholder="Enter your message"
                  variant="filled"
                  multiline
                  rows={8}
                  rowsMax={20}
                  value={contactDetails.message}
                  onChange={handleTextChange}
                  fullWidth
                  required
                  type="text"
                />
              </Grid>
              <Grid item>
                <Button
                  variant="contained"
                  style={{
                    height: '60px', width: '180px', borderRadius: '10px', fontWeight: 'bold',
                  }}
                  size="large"
                  onClick={formSubmit}
                  color="primary"
                >
                  Send
                </Button>
              </Grid>

            </Grid>
          </Box>
          <Footer />
        </div>
      </div>
    </ThemeProvider>

  );
};

export default Contact;
