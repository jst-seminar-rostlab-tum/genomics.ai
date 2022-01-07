import './app.module.css';
import React from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import HomePage from './components/LandingPage/pages/Home/Home';
import About from './components/LandingPage/pages/About/About';
import Docs from './components/LandingPage/pages/Docs/Docs';
import Contact from './components/LandingPage/pages/Contact/Contact';
import PasswordResetPage from './components/LandingPage/PasswordReset/PasswordResetPage';

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
          <Route exact path="/" component={HomePage} />
          <Route path="/about" component={About} />
          <Route path="/docs" component={Docs} />
          <Route path="/contact" component={Contact} />
          <Route path="/passwordreset" component={PasswordResetPage} />
        </Switch>
      </Router>
    </ThemeProvider>
  );
}

export default App;

// removed the dashboar from the router
// <Route path="/dashboard" component={Dashboard} />
// <Route exact path="/" component={LandingPage} />
