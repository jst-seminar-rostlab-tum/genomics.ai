import { Box } from '@mui/material';
import React, { useCallback } from 'react';
import {
  Redirect, Route, Switch, useRouteMatch,
} from 'react-router-dom';
import GeneMapperState from './GeneMapperState';
import GeneMapperHome from './Home';
import GeneMapperResultView from './Result';

function GeneMapper({ sidebarShown }) {
  const paddingL = useCallback(() => (sidebarShown ? '100px' : '350px'), [sidebarShown]);

  const { path, url } = useRouteMatch();

  return (
    <Box sx={{
      pl: paddingL,
      pr: '20px',
      height: '100vh',
      width: '100vw',
    }}
    >
      <Switch>
        <Route exact path={`${path}`}>
          <GeneMapperHome basePath={path} />
        </Route>

        <Route path={`${path}/create`}>
          <GeneMapperState path={`${path}`} />
        </Route>

        <Route path={`${path}/result/:projectId`}>
          <GeneMapperResultView />
        </Route>

        <Route path={`${path}/*`}>
          <Redirect to={`${url}`} />
        </Route>
      </Switch>
    </Box>
  );
}

export default GeneMapper;
