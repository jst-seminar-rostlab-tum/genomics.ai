import './index.css';
import { BrowserRouter as Router, Switch, Route } from 'react-router-dom';
import React from 'react';
import Stack from '@mui/material/Stack';
import Sidebar from './components/Dashboard/Sidebar/Sidebar';
import NavigationBar from './components/Dashboard/NavigationBar/NavigationBar';
import Dashboard from './components/Dashboard/Dashboard';
import Documentation from './components/Dashboard/Pages/Documentation/Documentation';
import Settings from './components/Dashboard/Pages/Settings/Settings';
import Help from './components/Dashboard/Pages/Help/Help';
import PopUp from './components/Dashboard/ProcessingCompletionPopUp/PopUp';

function App() {
  return (
    <Router>

      <NavigationBar />
      <Stack
        direction="row"
      >
        <Sidebar />
        <Switch>
          <Route path="/dashboard" exact component={Dashboard} />
          <Route path="/documentation" exact component={Documentation} />
          <Route path="/help" exact component={Help} />
          <Route path="/settings" exact component={Settings} />
        </Switch>
      </Stack>

      <PopUp/>
    </Router>
  );
}

export default App;

/*
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
*/
