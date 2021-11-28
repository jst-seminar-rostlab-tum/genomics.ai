import './app.module.css';
import React from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import LandingPage from './components/LandingPage/LandingPage';
import HomePage from './components/LandingPage/pages/HomePage/HomePage';
import About from './components/LandingPage/pages/About/About';
import Docs from './components/LandingPage/pages/Docs/Docs';

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
        </Switch>
      </Router>
    </ThemeProvider>
  );
}

export default App;

//removed the dashboar from the router
//<Route path="/dashboard" component={Dashboard} />
//<Route exact path="/" component={LandingPage} />