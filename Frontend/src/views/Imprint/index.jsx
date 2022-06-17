/* eslint-disable max-len */
import React, { useState, useCallback } from 'react';
import Box from '@mui/material/Box';
import Link from '@mui/material/Link';
import Typography from '@mui/material/Typography';
import { useHistory } from 'react-router-dom';
import NavBar from 'components/NavBar';
import Footer from 'components/Footer';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import { useAuth } from 'shared/context/authContext';

export default function Imprint() {
  const [isLoginFormVisible, setLoginFormVisible] = useState(false);
  const [isRegistrationFormVisible, setRegistrationFormVisible] = useState(false);

  const history = useHistory();
  const [user, setUser] = useAuth()

  const onLoginClicked = useCallback(() => {
    console.log('login');
    setRegistrationFormVisible(false);
    setLoginFormVisible(true);
  }, [setLoginFormVisible]);

  const onSignUpClicked = useCallback(() => {
    console.log('register');
    setLoginFormVisible(false);
    setRegistrationFormVisible(true);
  }, [setRegistrationFormVisible]);

  const onLoginFormClosed = useCallback(() => {
    setLoginFormVisible(false);
  }, [setLoginFormVisible]);

  const onRegistrationFormClosed = useCallback(() => {
    setRegistrationFormVisible(false);
  }, [setRegistrationFormVisible]);

  const executeScroll = () => user ? history.push({ pathname: '/sequencer/help'}) : history.push({ pathname: '/', state: { contact_us: true } });

  const regForm = isRegistrationFormVisible
    && <RegistrationForm visible={isRegistrationFormVisible} onClose={onRegistrationFormClosed} />;
  return (
    <Box sx={{ height: '100vh' }}>
      {isLoginFormVisible && <LoginForm visible={isLoginFormVisible} onClose={onLoginFormClosed} />}
      {regForm}

      <Box>
        <NavBar
          position="relative"
          onLoginClicked={onLoginClicked}
          onSignUpClicked={onSignUpClicked}
          executeScroll={executeScroll}
        />
      </Box>

      <Box sx={{
        margin: '3em auto',
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        width: {
          xs: '90%', sm: '90%', md: '70%', lg: '70%', xl: '70%',
        },
      }}
      >
        <Typography fontWeight="bold" fontSize="1.4em">Imprint</Typography>
      </Box>

      <Box sx={{
        margin: '3em auto',
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        width: {
          xs: '90%',
          sm: '90%',
          md: '70%',
          lg: '70%',
          xl: '70%',
        },
      }}
      >
        <Typography fontSize="1em">
          Helmholtz Zentrum München
          <br />
          Deutsches Forschungszentrum für Gesundheit und Umwelt (GmbH)
          <br />
          Ingolstädter Landstraße 1,
          <br />
          D-85764 Neuherberg, Germany
          <br />
          Phone: +49 89 3187-0
          <br />
          website: www.helmholtz-muenchen.de
          <br />
          email: presse(at)helmholtz-muenchen.de
          <br />
          <br />
          Representation-entitled board of directors: [placeholder]
          <br />
          Authorized representatives: [placeholder]
          <br />
          Register of Societies: Amtsgericht München HRB 6466
          <br />
          <br />
          VAT ID number in accordance with § 27 a Umsatzsteuergesetz (German Turnover-Tax Law): [placeholder]
          <br />
          <br />
          <br />
          <h1>Disclaimer</h1>
          <h2>1. Content</h2>

          The Helmholtz Zentrum München reserves the right not to be responsible for the topicality, correctness, completeness or quality of the information provided. Liability claims regarding damage caused by the use of any information provided, including any kind of information which is incomplete or incorrect will therefore be rejected.
          <br />
          Parts of the pages or the complete publication including all offers and information might be extended, changed or partly or completely deleted by the Helmholtz Zentrum München without separate announcement.

          <h2>2. Referrals and links</h2>

          The Helmholtz Zentrum München is not responsible for any contents linked or referred to from the Helmholtz Zentrum München´s pages - unless the Helmholtz Zentrum München has full knowledge of illegal contents and is able to prevent the visitors of the Helmholtz Zentrum München site from viewing those pages. If any damage occurs by the use of information presented there, only the author of the respective page might be liable, not the one who has linked to these pages.

          <h2>3. Copyright</h2>

          The Helmholtz Zentrum München intended not to use any copyrighted material for the publication or, if not possible, to indicate the copyright of the respective object.
          <br />
          The copyright for any material created by the Helmholtz Zentrum München is reserved. Any duplication or use of objects such as texts or diagrams, in other electronic or printed publications is not permitted without the Helmholtz Zentrum München´s prior agreement (please contact the public relation department: presse(at)helmholtz-muenchen.de).

          <h2>4. Privacy policy and data protection</h2>
          Our privacy policy can be found&nbsp;
          <Link href="/#/privacy">here</Link>
          .

          <h2>5. Legal validity of this disclaimer</h2>

          This disclaimer is to be regarded as part of the internet publication which you were referred from. If sections or individual terms of this statement are not legal or correct, the content or validity of the other parts remain unaffected by this fact.
        </Typography>
      </Box>

      <div style={{ height: 'calc(100vh - 485px)' }} />

      <Footer />
    </Box>
  );
}
