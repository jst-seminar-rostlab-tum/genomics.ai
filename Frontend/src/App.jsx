import './app.module.css';
import React, { useState } from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import {
  Route, HashRouter as Router, Switch, Redirect,
} from 'react-router-dom';
import HomePage from './views/Home';
import About from './views/About';
import Docs from './views/Docs';
import Contact from './views/Contact';
import DashboardContent from './components/DashboardContent';
import { guardedPage } from './shared/utils/common/utils';
import VisualizationPage from './views/VisualizationPage';
import PasswordResetPage from './views/PasswordResetPage';

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
      {/* eslint-disable-next-line no-restricted-globals */}
      <Router history={history}>
        <Switch>
          <Route exact path="/" render={() => (user ? <Redirect to="/sequencer" /> : <HomePage setUser={setUser} />)} />
          <Route path="/sequencer" render={() => guardedPage(<DashboardContent user={user} setUser={setUser} />)} />
          <Route path="/dashboard" component={DashboardContent} />
          <Route path="/about" render={() => <About setUser={setUser} />} />
          <Route path="/docs" render={() => <Docs setUser={setUser} />} />
          <Route path="/contact" render={() => <Contact setUser={setUser} />} />
          <Route path="/password_reset" render={() => <PasswordResetPage />} />
          <Route path="/result" render={() => <VisualizationPage />} />
        </Switch>
      </Router>
    </ThemeProvider>
  );
}

export default App;
