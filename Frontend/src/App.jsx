import './app.module.css';
import React from 'react';
import { ThemeProvider } from '@mui/material';
import {
  Route, HashRouter, Switch, Redirect,
} from 'react-router-dom';
import HomePage from './views/Home';
import About from './views/About';
import LegalNotice from './views/LegalNotice';
import Terms from './views/Terms';
import Privacy from './views/Privacy';
import Docs from './views/Docs';
import Contact from './views/Contact';
import DashboardContent from './components/DashboardContent';
import { guardedPage } from './shared/utils/common/utils';
import VisualizationPage from './views/VisualizationPage';
import PasswordResetPage from './views/PasswordResetPage';
import NotFound from './views/NotFound';
import { theme } from './shared/theme/theme';
import References from './views/References/index';
import { useAuth } from 'shared/context/authContext';
import NonLoginGeneMapper from 'views/NonLoginGeneMapper';

function App() {
  const [user] = useAuth();

  return (
    <ThemeProvider theme={theme}>
      {/* eslint-disable-next-line no-restricted-globals */}
      <HashRouter>
        <Switch>
          <Route exact path="/" render={() => (user ? <Redirect to="/sequencer" /> : <HomePage />)} />
          <Route path="/sequencer" render={() => guardedPage(<DashboardContent />)} />
          <Route path="/dashboard" component={DashboardContent} />
          {/* <Route path="/about" render={() => <About />} /> */}
          <Route path="/docs" render={() => <Docs />} />
          <Route path="/contact" render={() => <Contact />} />
          <Route path="/password_reset" render={() => <PasswordResetPage />} />
          <Route path="/result" render={() => <VisualizationPage />} />
          <Route path="/references" render={() => <References />} />
          {/* Gene mapper page without login */}
          <Route path="/genemapper" render={() => <NonLoginGeneMapper />} />
          <Route path="/imprint" render={() => <LegalNotice />} />
          <Route path="/terms" render={() => <Terms />} />
          <Route path="/privacy" render={() => <Privacy />} />
          <Route path="*" render={()=> <NotFound />}/>
        </Switch>
      </HashRouter>
    </ThemeProvider>
  );
}

export default App;
