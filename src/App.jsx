import './app.module.css';
import React from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import LandingPage from './components/LandingPage/LandingPage';
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
      <Router>
        <Switch>
          <Route exact path="/" component={LandingPage} />
          <Route path="/dashboard" component={Dashboard} />
        </Switch>
      </Router>
    </ThemeProvider>
  );
}

export default App;
