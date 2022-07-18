import { useContext, React } from 'react';
import { useHistory } from 'react-router-dom';
import { Box } from '@mui/material';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import PasswordForgetForm from 'components/PasswordForgetForm';
import NavBar from 'components/NavBar';
import Footer from 'components/Footer';
import { LoginContext } from 'shared/context/loginContext';
import GeneMapper from 'views/GeneMapper';
import AtlasModelChoice from 'views/GeneMapper/AtlasModelChoice/AtlasModelChoice';
import GeneMapperState from 'views/GeneMapper/GeneMapperState';

// The non-login version of the gene mapper page
function NonLoginGeneMapper() {
  const context = useContext(LoginContext);

  const onLoginClicked = () => {
    context.switchRegister(false);
    context.switchLogin(true);
  };

  const onSignUpClicked = () => {
    context.switchLogin(false);
    context.switchRegister(true);
  };

  const history = useHistory();

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', justifyContent: 'space-between' }}>
      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          minHeight: '100vh',
        }}
      >
        {context.loginVisible && <LoginForm />}
        {context.registerVisible && <RegistrationForm />}
        {context.forgetVisible && <PasswordForgetForm />}

        <Box>
          <NavBar
            position="relative"
            onLoginClicked={onLoginClicked}
            onSignUpClicked={onSignUpClicked}
            executeScroll={() => history.pushState({ pathname: '/', state: { contact_us: true } })} /* the history.push here might be wrong here */
          />
        </Box>
        <Box>
          {/* <GeneMapper sidebarShown={false} /> */}
          <AtlasModelChoice />
          {/* <GeneMapperState path={""}/> */}
        </Box>
      </Box>

      <Footer sx={{ marginTop: 'auto', transform: 'translate(0px, 100px)' }} />
    </Box>
  );
}

export default NonLoginGeneMapper;
