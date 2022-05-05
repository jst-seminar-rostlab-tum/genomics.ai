import React, { useCallback, useState } from 'react';
import {
  Stack, ThemeProvider, createTheme, Typography, TextField, Button,
} from '@mui/material';
import styles from './help.module.css';

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

function Documentation() {
  const [contactDetails, setContactDetails] = useState({
    message: '',
  });

  const handleTextChange = useCallback((e) => {
    setContactDetails((prevState) => ({ ...prevState, [e.target.id]: e.target.value }));
  }, [setContactDetails]);

  const formSubmit = async (e) => {
    e.preventDefault();

    const data = {
      message: contactDetails.message,
    };
    console.log(data);

    try {
      // send POST to the server
      // await axios.post("BACKEND_URL", data);
      setContactDetails({
        message: '...sending',
      });
    } catch (error) {
      console.log(error);
    }
  };

  return (

    <ThemeProvider theme={theme}>
      <Stack
        direction="column"
        sx={{
          paddingTop: '100px',
          paddingLeft: '130px',
          paddingRight: '30px',
          display: 'flex',
          width: '100%',
          justifyContent: 'center',
        }}
        spacing={1}
      >
        <div className={styles.title}>
          <h1>Help</h1>
        </div>

        <Stack direction="column" spacing={4}>
          <Stack spacing={0.5}>
            <Typography sx={{ fontSize: '20px' }}>
              Please describe to us the problem you are experiencing.
            </Typography>
            <Typography sx={{ fontSize: '20px' }}>
              We will get back to you within the next 48 hours.
            </Typography>
          </Stack>

          <TextField
            id="message"
            label="Please describe your problem"
            placeholder="Enter your message"
            variant="filled"
            multiline
            rows={8}
            rowsMax={20}
            value={contactDetails.message}
            onChange={handleTextChange}
            required
            type="text"
            sx={{ maxWidth: '1000px' }}
          />

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
        </Stack>

      </Stack>

    </ThemeProvider>

  );
}

export default Documentation;
