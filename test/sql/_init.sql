\set ECHO none

-- disable the notices for the create script (shell types etc.)
SET client_min_messages = 'WARNING';
\i sql/ddsketch--1.0.0.sql
\i sql/ddsketch--1.0.0--1.0.1.sql
\i sql/ddsketch--1.0.1--1.0.2-dev.sql
SET client_min_messages = 'NOTICE';
