import './app.module.css';
import React, { useState } from 'react';
import { ThemeProvider } from '@mui/material';
import {
  Route, HashRouter, Switch, Redirect,
} from 'react-router-dom';
import HomePage from './views/LandingPage/Home';
import About from './views/LandingPage/About';
import Docs from './views/LandingPage/Docs';
import Contact from './views/LandingPage/Contact';
import DashboardContent from './components/DashboardContent';
import { guardedPage } from './shared/utils/common/utils';
import VisualizationPage from './views/VisualizationPage';
import PasswordResetPage from './views/LandingPage/PasswordResetPage';
import { theme } from "./shared/theme/theme"
import Explore from "./views/Explore/index.jsx"

function App() {

  const [user] = useAuth()

  return (
    <ThemeProvider theme={theme}>
      {/* eslint-disable-next-line no-restricted-globals */}
      <HashRouter>
        <Switch>
          <Route exact path="/" render={() => (user ? <Redirect to="/sequencer" /> : <HomePage />)} />
          <Route path="/sequencer" render={() => guardedPage(<DashboardContent />)} />
          <Route path="/dashboard" component={DashboardContent} />
          <Route path="/about" render={() => <About />} />
          <Route path="/docs" render={() => <Docs />} />
          <Route path="/contact" render={() => <Contact />} />
          <Route path="/password_reset" render={() => <PasswordResetPage />} />
          <Route path="/result" render={() => <VisualizationPage />} />
          <Route path="/explore" render={() => <Explore />} />
        </Switch>
      </HashRouter>
    </ThemeProvider>
  );
}

export default App;
