import { useContext, React } from 'react';
import { useHistory } from 'react-router-dom';
import { Box } from '@mui/material';
import GeneMapper from 'views/GeneMapper';
import LoginForm from 'components/LoginForm';
import RegistrationForm from 'components/RegistrationForm';
import PasswordForgetForm from 'components/PasswordForgetForm';
import { LoginContext } from 'shared/context/loginContext';

// The non-login version of the genemapper page.
// The main genemapper component is used with loggedIn set to false.
function NonLoginGeneMapper() {
  return (
    <Box sx={{
      display: 'flex',
      flexDirection: 'column',
      justifyContent: 'space-between',
      mineight: '100vh',
    }}
    >
      <Box
        sx={{
          display: 'flex',
          flexDirection: 'column',
          minHeight: '100vh',
        }}
      >
        {/* Non-login GeneMapper */}
        <GeneMapper sidebarShown={false} loggedIn={false} />
      </Box>
    </Box>
  );
}

export default NonLoginGeneMapper;
