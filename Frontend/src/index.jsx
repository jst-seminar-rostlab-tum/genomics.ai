import React from 'react';
import ReactDOM from 'react-dom';
import './index.css';
import { CssBaseline } from '@mui/material';
import App from './App';
import { AuthProvider } from "shared/context/authContext"

ReactDOM.render(
  <>
    <CssBaseline />
    <AuthProvider><App /> </AuthProvider>
  </>,
  document.getElementById('root'),
);
