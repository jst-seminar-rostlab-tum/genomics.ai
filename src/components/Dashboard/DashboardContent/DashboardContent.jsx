import React, { useState } from 'react';
import { BrowserRouter as Router, Switch, Route } from 'react-router-dom';
import Sidebar from '../Sidebar/Sidebar';
import Dashboard from '../Dashboard';
import NavigationBar from '../NavigationBar/NavigationBar';
import Documentation from '../Pages/Documentation/Documentation';
import Settings from '../Pages/Settings/Settings';
import Help from '../Pages/Help/Help';
import PopUp from '../ProcessingCompletionPopUp/PopUp';

const DashboardContent = () => {
  const [sidebarShown, setSidebarShown] = useState(true);
  const toggleSidebar = () => setSidebarShown(!sidebarShown);

  return (
    <Router>
      <NavigationBar sidebarShown={sidebarShown} />
      <Sidebar
        toggleSidebar={toggleSidebar}
        sidebarShown={sidebarShown}
      />

      <Switch>
        <Route path="/dashboard">
          <Dashboard sidebarShown={sidebarShown} />
        </Route>

        <Route path="/documentation">
          <Documentation sidebarShown={sidebarShown} />
        </Route>

        <Route path="/help">
          <Help sidebarShown={sidebarShown} />
        </Route>

        <Route path="/settings">
          <Settings sidebarShown={sidebarShown} />
        </Route>

      </Switch>

      <PopUp />

    </Router>

  );
};

export default DashboardContent;
