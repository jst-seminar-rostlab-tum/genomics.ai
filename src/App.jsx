import './app.module.css';
import React from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import LandingPage from './components/LandingPage/LandingPage';

function App() {
  const theme = createTheme({
    palette: {
      primary: {
        main: '#01579B',
      },
      light: {
        main: '#4F83CC',
      },
    },
  });

  return (
    <ThemeProvider theme={theme}>
      <div>
        <LandingPage />
      </div>
    </ThemeProvider>
  );
}

export default App;
