import './app.module.css';
import React, { useState } from 'react';
import { ThemeProvider, createTheme } from '@mui/material';
import {
  Route, HashRouter, Switch, Redirect,
} from 'react-router-dom';
import HomePage from './views/Home';
import About from './views/About';
import Docs from './views/Docs';
import Contact from './views/Contact';
import DashboardContent from './components/DashboardContent';
import { guardedPage } from './shared/utils/common/utils';
import VisualizationPage from './views/VisualizationPage';
import PasswordResetPage from './views/PasswordResetPage';
import { theme } from "./shared/theme/theme"
import Explore from "./views/Explore/index.jsx"
import UploadFilePage from 'views/GeneMapper/UploadFilePage';

function App() {
  // https://stackoverflow.com/a/69836010
  const { palette } = createTheme();
  const { augmentColor } = palette;
  const createColor = (mainColor) => augmentColor({ color: { main: mainColor } });
  const theme = createTheme({
    palette: {
      primary: {
        main: '#01579B',
      },
      light: {
        main: '#4F83CC',
      },
      critical: createColor('#F44336'),
    },
  });

  const [user, setUser] = useState(localStorage.user ? JSON.parse(localStorage.user) : null);

  return (
    <ThemeProvider theme={theme}>
      {/* eslint-disable-next-line no-restricted-globals */}
      <HashRouter>
        <Switch>
          <Route exact path="/" render={() => (user ? <Redirect to="/sequencer" /> : <HomePage setUser={setUser} />)} />
          <Route path="/sequencer" render={() => guardedPage(<DashboardContent user={user} setUser={setUser} />)} />
          <Route path="/dashboard" component={DashboardContent} />
          <Route path="/about" render={() => <About setUser={setUser} />} />
          <Route path="/docs" render={() => <Docs setUser={setUser} />} />
          <Route path="/contact" render={() => <Contact setUser={setUser} />} />
          <Route path="/password_reset" render={() => <PasswordResetPage />} />
          <Route path="/result" render={() => <VisualizationPage />} />
          <Route path="/explore" render={() => <Explore />} />
          <Route path="/genemapper" render={() => <UploadFilePage />} />
        </Switch>
      </HashRouter>
    </ThemeProvider>
  );
}

export default App;
