import React, { useState, useCallback } from 'react';
import {
  Typography, Divider, Stack, createTheme, ThemeProvider, Button,
} from '@mui/material';
import NavBar from '../../../NavBar/NavBar';
import styles from './home.module.css';
import graphic1 from '../../../../assets/landing-illustrations/science.png';
import graphic2 from '../../../../assets/landing-illustrations/upload.png';
import graphic3 from '../../../../assets/landing-illustrations/processing.png';
import graphic4 from '../../../../assets/landing-illustrations/results.png';
import tum from '../../../../assets/landing-illustrations/tum-logo.png';
import rostlab from '../../../../assets/landing-illustrations/rostlab.png';
import helmholtz from '../../../../assets/landing-illustrations/helmholtz.png';
import Footer from '../../Footer/Footer';
import dnaImage from '../../../../assets/dna.png';
import RegistrationForm from '../../RegistrationForm/RegistrationForm';
import LoginForm from '../../LoginForm/LoginForm';

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

function Home() {
// inserting the logic from the landing page here
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const onLoginClicked = useCallback(() => {
    setRegistrationFormVisible(false);
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

  const onSignUpClicked = useCallback(() => {
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible]);

  const onLoginFormClosed = useCallback(() => {
    setLoginFormVisible(false);
  }, [setLoginFormVisible]);

  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  return (
    <div className={styles.container}>
      <NavBar />
      {/* If genomics visualized stays in the middle, then it needs to be properly aligned.
      Right now it is using padding as the way to align the items */}

      <div className={styles.call2Action}>
        <Typography
          sx={{
            fontSize: '55px', fontWeight: 'bold', color: 'white', paddingInline: '150px',
          }}
        >
          AI-driven Cell Type Annotation
        </Typography>
        <Typography
          sx={{
            fontSize: '55px', fontWeight: 'regular', color: 'white', paddingInline: '150px',
          }}
        >
          No Code. Just Results.
        </Typography>
        <Typography
          sx={{
            fontSize: '25px', fontWeight: 'regular', color: 'white', paddingBlock: '40px', paddingInline: '300px',
          }}
        >
          Genomics.ai helps you visualize all of your single-cell sequencing data
          in a fast and easy way using neural networks.
        </Typography>

        <ThemeProvider theme={theme}>
          <Stack
            spacing={2}
            direction="row"
            style={{ height: '60px', marginTop: '30px', justifyContent: 'center' }}
            size="large"
          >
            <Button
              variant="contained"
              color="primary"
              style={{
                width: '200px', borderRadius: '12px', fontWeight: 'bold',
              }}
              size="large"
              onClick
              onClick={onSignUpClicked}
            >
              Sign Up (FREE)
            </Button>

            <Button
              variant="outlined"
              color="secondary"
              style={{ width: '150px', borderRadius: '12px' }}
              href="contact"
            >
              Get in Touch
            </Button>

            <div>
              <LoginForm
                visible={isLoginFormVisible}
                onClose={onLoginFormClosed}
              />
              <RegistrationForm
                visible={isRegistrationFormVisible}
                onClose={onRegistrationFormClosed}
                onSuccessfulRegistration={onLoginClicked}
              />

            </div>
          </Stack>
        </ThemeProvider>
      </div>

      <img
        className={styles.backgroundImage}
        src={dnaImage}
        alt="picture of DNA"
      />

      <Divider
        variant="middle"
        textAlign="center"
        sx={{ padding: '30px', paddingBottom: '50px', paddingTop: '100px' }}
      >
        <Typography sx={{ fontSize: '35px', fontWeight: 'bold' }}>Our Partners</Typography>
      </Divider>

      <Stack
        container
        direction="row"
        justifyContent="space-evenly"
        alignItems="center"
        padding="200px"
        paddingTop="80px"
      >
        <img className={styles.tumLogo} src={tum} alt="tum-logo" />
        <img className={styles.rostlabLogo} src={rostlab} alt="rostlab" />
        <img className={styles.tumLogo} src={helmholtz} alt="helmholtz" />
      </Stack>

      <Divider
        variant="middle"
        textAlign="center"
        sx={{ padding: '30px', paddingTop: '50px' }}
      >
        <Typography sx={{ fontSize: '35px', fontWeight: 'bold' }}>What we do</Typography>
      </Divider>

      <div className={styles.infoContainer}>
        <img className={styles.illustration} src={graphic1} alt="science-guy" />
        <div className={styles.explanation}>
          <Typography
            className={styles.illustrationTitle}
            sx={{ fontSize: '30px', fontWeight: 'bold' }}
          >
            Genomics.ai
          </Typography>
          <Typography sx={{ fontSize: '25px' }}>
            {' '}
            We help you visualize all
            of your single-cell sequencing data in a fast and easy way with the help of neural networks.
          </Typography>
        </div>
      </div>

      <Divider variant="middle" textAlign="center" sx={{ padding: '30px', paddingBottom: '50px', paddingTop: '50px' }}>
        <Typography sx={{ fontSize: '35px', fontWeight: 'bold' }}>How it works</Typography>
      </Divider>

      <div className={styles.infoContainer}>
        <div className={styles.explanation}>
          <Typography className={styles.illustrationTitle} sx={{ fontSize: '30px', fontWeight: 'bold' }}>Upload</Typography>
          <Typography sx={{ fontSize: '25px' }}>Upload your files containing the single-cell sequencing data.</Typography>
        </div>
        <img className={styles.illustration} src={graphic2} alt="upload" />
      </div>

      <div className={styles.infoContainer}>
        <img className={styles.illustration} src={graphic3} alt="processing" />
        <div className={styles.explanation}>
          <Typography
            className={styles.illustrationTitle}
            sx={{ fontSize: '30px', fontWeight: 'bold' }}
          >
            Processing
          </Typography>
          <Typography sx={{ fontSize: '25px' }}>
            Your input data is processed by our machine learning model
            which performs the cell-type classification according to your specification.
          </Typography>
        </div>
      </div>

      <div className={styles.infoContainer}>
        <div className={styles.explanation}>
          <Typography
            className={styles.illustrationTitle}
            sx={{ fontSize: '30px', fontWeight: 'bold' }}
          >
            Check Results
          </Typography>
          <Typography sx={{ fontSize: '25px' }}>
            After the algorithm has processed the data you can
            specify the project you want your cell-type data to be associated with and view the results.
          </Typography>
        </div>
        <img className={styles.illustration} src={graphic4} alt="results" />
      </div>

      <Footer />
    </div>
    // TODO: add footer for the website across all pages that are not the dashboard

  );
}

export default Home;
