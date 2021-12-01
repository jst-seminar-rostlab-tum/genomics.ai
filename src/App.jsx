import './app.module.css';
import React from 'react';
import { createTheme, ThemeProvider } from '@mui/material';
import { Route, BrowserRouter as Router, Switch } from 'react-router-dom';
import HomePage from './components/LandingPage/pages/Home/Home.jsx';
import About from './components/LandingPage/pages/About/About';
import Docs from './components/LandingPage/pages/Docs/Docs';
import Contact from './components/LandingPage/pages/Contact/Contact';
import Dashboard from './components/Dashboard/Dashboard';
import DashboardContent from './components/Dashboard/DashboardContent/DashboardContent';
import Visualization from './components/Dashboard/Visualization/src/Visualization';

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
          <Route path="/dashboard" component={DashboardContent} />
          <Route path="/about" component={About} />
          <Route path="/docs" component={Docs} />
          <Route path="/contact" component={Contact} />
          <Route path="/visualization" component={Visualization} />
        </Switch>
      </Router>
    </ThemeProvider>
  );
}

export default App;

// removed the dashboar from the router
// <Route path="/dashboard" component={Dashboard} />
// <Route exact path="/" component={LandingPage} />
