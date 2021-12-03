import './app.module.css';
import React, { useState } from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import {
  Route, HashRouter as Router, Switch, Redirect,
} from 'react-router-dom';
import HomePage from './components/LandingPage/pages/Home/Home';
import About from './components/LandingPage/pages/About/About';
import Docs from './components/LandingPage/pages/Docs/Docs';
import Contact from './components/LandingPage/pages/Contact/Contact';
import DashboardContent from './components/Dashboard/DashboardContent/DashboardContent';
import { guardedPage } from './components/common/utils';

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

  const [user, setUser] = useState(localStorage.user ? JSON.parse(localStorage.user) : null);

  return (
    <ThemeProvider theme={theme}>
      <Router>
        <Switch>
          <Route exact path="/" render={() => (user ? <Redirect to="/sequencer" /> : <HomePage setUser={setUser} />)} />
          <Route path="/sequencer" render={() => guardedPage(<DashboardContent user={user} setUser={setUser} />)} />
          <Route path="/dashboard" component={DashboardContent} />
          <Route path="/about" render={() => <About setUser={setUser} />} />
          <Route path="/docs" render={() => <Docs setUser={setUser} />} />
          <Route path="/contact" render={() => <Contact setUser={setUser} />} />
        </Switch>
      </Router>
    </ThemeProvider>
  );
}

export default App;
