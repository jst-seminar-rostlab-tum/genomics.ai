import './app.module.css';
import React from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import Dashboard from './components/Dashboard/Dashboard';

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
        <Dashboard />
      </div>
    </ThemeProvider>
  );
}

export default App;
